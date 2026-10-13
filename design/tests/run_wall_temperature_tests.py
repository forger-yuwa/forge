#!/usr/bin/env python3
"""壁温分布 $T_w(x)$ の配管 (plan tooling-nozzle-wall-temperature-profile §5.1 #2) の試験 — §6 V1・V1b・V1d・V4。

  V1  書式の検査 (`wall_temperature.parse_profile` の正例 s・x_m と負例: 非単調・負・NaN・Inf・長さ違い・座標キー 2 つ・
      s が [−2, 1] を覆わない・`Tw` との併用 ほか)、区間の写像 (`to_x_m`) の数値が手計算と相対 1e-15、`Tw_table` の x が
      x_m / r_t と相対 1e-15、bcond の壁の行に `wallProfile: 1`・`profile_ok=False` で例外、CSV の往復 (相対 1e-12)。
  V1b 実効の r_t の A/B (codex 諮問 2026-10-13): 問題と表を固定し `integral_delta_r` の scale だけを A = spec.r_throat・B = 1.2 A。
      積分法が実際に使った温度 res["Tw"] と、原表を実寸座標 scale × res["x"] で独立に評価した温度の最大差が A・B とも ≤ 1e-8 K。
      x_m の表で B だけ積分区間が表の外に出ると例外。旧半径で換算する誤実装なら B で 10 K 以上ずれる (表の勾配の確認)。
      寸法の逆算の CFD 前の供給 (`deltastar_loop._sizing_delta_supplier`) も同じ scale で表を換算する。
  V1d 正本の照合: `deltastar_loop --init-integral --prepare-only` で作った run の CSV を 1 値書き換える → `load_run_profile` が止まり、
      `run_staged_ns` が forge を起動する前に止まる (FORGE_BIN を存在しない道にして、起動の痕跡 [run_case_stdout.log・
      RUN_PROVENANCE.txt・forge_run.log] が無いことで確かめる)。スロートの温度が同じで分布だけ違う表で作った別の run の JSON を
      写すと止まる。段の間 (restart の後) に CSV を書き換えても次の段の起動前に止まる。書き換えていない run は照合を通り起動まで進む (対照)。
  V4  後方互換: case/45 の生産 YAML (断熱) と case/70 の c2b (一様 350 K) で、`_bcond(p, euler=False)` の文字列・
      `integral_thermal_bc` の dict が変更前 (REF_COMMIT) の `wall_bcond_line`・`wall_thermal_bc_integral` の結果と同一。
      さらに変更前のコード (REF_COMMIT の design/・solver_density_cuda/data を git archive) と今のコードで同じ問題を `--prepare-only` し、生成物
      (bcondConfig・solverConfig・delta_r_initial ほか・prepare_info から environment と壁の当てはめの所要時間を除いたもの) がバイト同一。
      今のコードの側は逆 MOC の出口の切り方を直す前の規則に戻して回す (2026-10-13 plan tooling-nozzle-deltastar-band-root-zone §5.1 #4c で
      意図して変えた; 両問題とも直す前の切り方が x_F の壁点を落としていた)。今の切り方では壁の長さに依るファイルだけが違い、設計壁の終点が
      下流へ動くことを別に見る。

forge の計算は回さない (変換器と `forge --resolve-species` は prepare が使う)。出力は WALL_TEMP_TEST_OUT (無ければ一時ディレクトリで、
全部合格なら最後に消す; run の写しで 700 MB ほどになる)。
既存の run (case/45・case/70) は読むだけ。
実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_wall_temperature_tests.py
所要: 約 2 分 (prepare 6 回・積分法 6 回)。
"""
import ast
import copy
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
REPO = DESIGN.parent
BUILD = REPO / "solver_density_cuda" / ".build-native" / "relwithdebinfo"
# prepare が使う変換器と forge (--resolve-species)。runner の _ENV は import 時の環境を写すので import の前に決める
os.environ.setdefault("FORGE_CONVERTER", str(BUILD / "convertGmshToForge"))
os.environ.setdefault("FORGE_BIN", str(BUILD / "forge"))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")
sys.path.insert(0, str(DESIGN))

from forge_design import wall_temperature as WTP  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate.runner_wt import _bcond  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

C70 = REPO / "case" / "70.m5_air_wt_demo" / "problem_m5_air_c2b.yaml"           # 一様 350 K
C70_RUN = REPO / "case" / "70.m5_air_wt_demo" / "run_0010_ns_c2b_final"           # 一様の既存 run (読むだけ)
C45 = REPO / "case" / "45.isobutane_m6_d155" / "problem_d155_ns_prod.yaml"        # 断熱 (生産)
C45_PIN_DIR = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")              # initial_line_run (CFD ピン) の所在 (読むだけ)
REF_COMMIT = "138e83b2dab18451eaddc35cd185075932d8f6f7"   # 壁温分布の配管を入れる前の commit (V4 の「変更前」)
KEEP_OUT = bool(os.environ.get("WALL_TEMP_TEST_OUT"))
OUT = Path(os.environ.get("WALL_TEMP_TEST_OUT") or tempfile.mkdtemp(prefix="wall_temp_tests_"))
OUT.mkdir(parents=True, exist_ok=True)
PY = sys.executable
NO_FORGE = "/nonexistent/forge_wall_temperature_test"
FAIL = 0
SUMMARY = {}

# 表 (内部で 10 K 以上の勾配)。s: 区間で正規化 / x_m: 実寸 [m] (case/70: 入口端 −0.500 m・壁の終点 1.1995 m、B では 1.2 倍)
PROF_S = {"s": [-2.0, -1.0, -0.5, -0.2, 0.0, 0.1, 0.4, 1.0], "Tw_K": [320.0, 340.0, 400.0, 470.0, 500.0, 480.0, 380.0, 330.0]}
PROF_S2 = {"s": [-2.0, -1.0, -0.3, 0.0, 0.3, 1.0], "Tw_K": [300.0, 310.0, 450.0, 500.0, 420.0, 350.0]}   # スロート 500 K は同じ
PROF_X = {"x_m": [-0.65, -0.3, -0.1, 0.0, 0.05, 0.3, 0.8, 1.5], "Tw_K": [320.0, 360.0, 450.0, 500.0, 480.0, 400.0, 340.0, 330.0]}
PROF_XN = {"x_m": [-0.52, -0.3, -0.1, 0.0, 0.05, 0.3, 0.8, 1.25], "Tw_K": [320.0, 360.0, 450.0, 500.0, 480.0, 400.0, 340.0, 330.0]}


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def raises(fn, *substr, exc=ValueError):
    """fn() が exc を出し、メッセージが substr を全部含むか。(結果, メッセージ)。"""
    try:
        fn()
    except exc as ex:
        msg = str(ex)
        return all(s in msg for s in substr), msg[:220]
    except Exception as ex:          # noqa: BLE001 (別の例外は不合格として内容を出す)
        return False, f"別の例外 {type(ex).__name__}: {str(ex)[:200]}"
    return False, "例外が出ない"


def rel_ok(a, b, tol):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    return a.shape == b.shape and bool(np.all(np.abs(a - b) <= tol * np.maximum(np.abs(b), 1e-300)))


def write_problem(base: Path, wall_thermal: dict, name: str, abs_pin: bool = False) -> Path:
    """base の問題 YAML を写し、spec.wall_thermal を置き換えて OUT に書く (abs_pin: initial_line_run を絶対パスに)。"""
    raw = yaml.safe_load(base.read_text())
    raw["spec"]["wall_thermal"] = copy.deepcopy(wall_thermal)
    if abs_pin and raw["geometry"].get("initial_line_run"):
        raw["geometry"]["initial_line_run"] = str(C45_PIN_DIR / raw["geometry"]["initial_line_run"])
    path = OUT / f"{name}.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    return path


def ref_x_of_table(prof: dict, scale: float, d: dict, geom: dict) -> np.ndarray:
    """試験側で独立に書いた写像: 表の点 → 実寸の x [m] (区切りは L_pipe・L_U・壁の終点 × scale)。"""
    if "x_m" in prof:
        return np.asarray(prof["x_m"], dtype=float)
    L_pipe, L_U = float(geom.get("L_pipe", 0.5)), float(geom.get("L_U", 3.5))
    xi, xc, xe = -(L_pipe + L_U) * scale, -L_U * scale, float(d["wall"].x_e) * scale
    out = []
    for s in prof["s"]:
        if s < -1.0:
            out.append(xi + (s + 2.0) * (xc - xi))
        elif s < 0.0:
            out.append(xc + (s + 1.0) * (0.0 - xc))
        else:
            out.append(s * xe)
    return np.asarray(out)


def prepare_cli(design_dir: Path, problem: Path, run_dir: Path, env_extra=None, old_exit_cut: bool = False):
    """`deltastar_loop --init-integral --prepare-only` (CFD 前の準備) を design_dir のコードで回す。戻り (rc, ログ)。
    `old_exit_cut`: 逆 MOC の出口の切り方を直す前 (plan tooling-nozzle-deltastar-band-root-zone §5.1 #4c) の規則 `wall[:, 0] <= x_F` に
    戻して回す (変更前のコードとの比較で、意図して変えた切り方の差を除くため)。"""
    if run_dir.exists():
        shutil.rmtree(run_dir)
    env = dict(os.environ, FORGE_ROOT=str(REPO), FORGE_NO_DESIGN_REPORT="1", **(env_extra or {}))
    args = ["--problem", str(problem), "--euler-ref", "none", "--run-dir", str(run_dir), "--init-integral", "--prepare-only"]
    if old_exit_cut:
        code = ("import runpy, sys; from forge_design.geometry import moc_inverse as MI; "
                "MI._cut_wall_at_exit = lambda w, xF: (w[w[:, 0] <= xF], {}); "
                f"sys.argv = ['deltastar_loop'] + {args!r}; "
                "runpy.run_module('forge_design.feedback.deltastar_loop', run_name='__main__')")
        cmd = [PY, "-c", code]
    else:
        cmd = [PY, "-m", "forge_design.feedback.deltastar_loop", *args]
    r = subprocess.run(cmd, cwd=design_dir, env=env, capture_output=True, text=True)
    (run_dir.parent / f"{run_dir.name}.prepare.log").write_text(r.stdout + r.stderr)
    return r.returncode, r.stdout + r.stderr


def staged_cli(run_dir: Path, stages: str = "full"):
    """`run_staged_ns` を FORGE_BIN = 存在しない道で子プロセスで呼ぶ。戻り (rc, ログ)。"""
    env = dict(os.environ, FORGE_BIN=NO_FORGE, FORGE_ROOT=str(REPO))
    code = ("import sys; from forge_design.evaluate.runner_axismach import run_staged_ns; "
            f"rc = run_staged_ns({str(run_dir)!r}, stages={stages!r}); print('RUN_STAGED_RC', rc)")
    r = subprocess.run([PY, "-c", code], cwd=DESIGN, env=env, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


LAUNCH_TRACES = ("run_case_stdout.log", "RUN_PROVENANCE.txt", "forge_run.log", "stage_manifest.json")


def copy_run(src: Path, name: str) -> Path:
    dst = OUT / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return dst


# ======================================================================================================
print(f"出力: {OUT}", flush=True)
P_S = load_problem(write_problem(C70, {"mode": "isothermal", "Tw_profile": PROF_S}, "c70_prof_s"))
P_S2 = load_problem(write_problem(C70, {"mode": "isothermal", "Tw_profile": PROF_S2}, "c70_prof_s2"))
P_X = load_problem(write_problem(C70, {"mode": "isothermal", "Tw_profile": PROF_X}, "c70_prof_x"))
P_XN = load_problem(write_problem(C70, {"mode": "isothermal", "Tw_profile": PROF_XN}, "c70_prof_xn"))
P_U = load_problem(C70)
D = RA.design_chain(P_U)                  # 設計は wall_thermal に依らない (同じ geometry なので全部の表で共有)
A = float(P_U.spec["r_throat"])
B = 1.2 * A
GEOM = P_U.geometry
print(f"case/70 c2b: r_t = {A!r} m、壁 x_in {D['wall'].x_in!r}・x_e {D['wall'].x_e!r} r_t", flush=True)

print("== V1 書式の検査 (parse_profile)", flush=True)
ps = WTP.parse_profile(PROF_S)
check("s の表 (正例)", ps == {"coord": "s", "coords": PROF_S["s"], "Tw_K": PROF_S["Tw_K"]}, str(ps)[:120])
px = WTP.parse_profile(PROF_X)
check("x_m の表 (正例)", px == {"coord": "x_m", "coords": PROF_X["x_m"], "Tw_K": PROF_X["Tw_K"]})
check("YAML から読んだ表 (整数を含む) も float の列", WTP.parse_profile(yaml.safe_load("{s: [-2, 0, 1], Tw_K: [300, 500, 330]}"))
      == {"coord": "s", "coords": [-2.0, 0.0, 1.0], "Tw_K": [300.0, 500.0, 330.0]})
NEG = [
    ("非単調", {"s": [-2.0, 0.0, -1.0, 1.0], "Tw_K": [300, 400, 500, 300]}, "狭義単調増加"),
    ("座標の重複", {"s": [-2.0, 0.0, 0.0, 1.0], "Tw_K": [300, 400, 500, 300]}, "狭義単調増加"),
    ("負の温度", {"s": [-2.0, 0.0, 1.0], "Tw_K": [300, -5.0, 300]}, "正の温度"),
    ("0 K", {"s": [-2.0, 1.0], "Tw_K": [0.0, 300]}, "正の温度"),
    ("座標に NaN", {"s": [-2.0, float("nan"), 1.0], "Tw_K": [300, 400, 300]}, "有限でない"),
    ("温度に NaN (YAML の .nan)", yaml.safe_load("{s: [-2, 0, 1], Tw_K: [300, .nan, 300]}"), "有限でない"),
    ("座標に Inf", {"x_m": [-1.0, 0.0, float("inf")], "Tw_K": [300, 400, 300]}, "有限でない"),
    ("温度に Inf", {"x_m": [-1.0, 0.0, 1.0], "Tw_K": [300, float("inf"), 300]}, "有限でない"),
    ("長さ違い", {"s": [-2.0, 0.0, 1.0], "Tw_K": [300, 400]}, "長さが違う"),
    ("座標キー 2 つ", {"s": [-2.0, 1.0], "x_m": [-1.0, 1.0], "Tw_K": [300, 300]}, "どちらか 1 つ"),
    ("座標キーなし", {"Tw_K": [300, 300]}, "どちらか 1 つ"),
    ("Tw_K なし", {"s": [-2.0, 1.0]}, "Tw_K"),
    ("s が −2 を覆わない", {"s": [-1.9, 0.0, 1.0], "Tw_K": [300, 400, 300]}, "[−2, 1]"),
    ("s が 1 を覆わない", {"s": [-2.0, 0.0, 0.99], "Tw_K": [300, 400, 300]}, "[−2, 1]"),
    ("x_m が設計スロート x = 0 を含まない", {"x_m": [0.1, 1.0], "Tw_K": [300, 300]}, "x = 0"),
    ("1 点", {"s": [0.0], "Tw_K": [300]}, "2 点以上"),
    ("文字列の数", {"s": ["-2", "1"], "Tw_K": [300, 300]}, "数値"),
    ("2 次元 (入れ子)", {"s": [[-2.0, 1.0]], "Tw_K": [[300, 300]]}, "数値"),
    ("真偽値", {"s": [-2.0, True], "Tw_K": [300, 300]}, "数値"),
    ("スカラー", {"s": 0.0, "Tw_K": 300}, "数値の列"),
    ("未知のキー", {"s": [-2.0, 1.0], "Tw_K": [300, 300], "unit": "K"}, "未知のキー"),
    ("mapping でない", [[-2.0, 300.0], [1.0, 300.0]], "mapping"),
]
for name, obj, sub in NEG:
    ok, msg = raises(lambda o=obj: WTP.parse_profile(o), sub)
    check(f"負例: {name} → ValueError", ok, msg)

print("== V1 問題 YAML (probdef.wall_thermal)", flush=True)
wt_s = P_S.wall_thermal
check("分布: mode isothermal・代表値 Tw = s = 0 の値 (500 K)・Tw_profile = parse_profile",
      wt_s == {"mode": "isothermal", "Tw": 500.0, "Tw_profile": ps}, str(wt_s)[:160])
check("分布 (x_m): 代表値 = x = 0 の値 (500 K)", P_X.wall_thermal["Tw"] == 500.0)
for name, wtd, sub in (
        ("Tw と Tw_profile の併用", {"mode": "isothermal", "Tw": 350.0, "Tw_profile": PROF_S}, "排他"),
        ("Tw_profile_csv (カード専用) を問題 YAML に", {"mode": "isothermal", "Tw_profile_csv": "t.csv"}, "カード"),
        ("断熱に Tw_profile", {"mode": "adiabatic", "Tw_profile": PROF_S}, "isothermal"),
        ("壊れた表 (非単調)", {"mode": "isothermal", "Tw_profile": {"s": [-2.0, 1.0, 0.5], "Tw_K": [300, 300, 300]}}, "狭義単調増加")):
    pn = load_problem(write_problem(C70, wtd, "neg_" + str(abs(hash(name)) % 10 ** 8)))
    ok, msg = raises(lambda q=pn: q.wall_thermal, sub)
    check(f"負例 (問題 YAML): {name} → ValueError", ok, msg)
ok, msg = raises(lambda: P_S.wall_thermal_bc_integral, "integral_thermal_bc")
check("分布で wall_thermal_bc_integral (scale を知らない property) は例外", ok, msg)

print("== V1 区間の写像 (to_x_m)・範囲の検査・補間", flush=True)
anc = WTP.anchors_m(0.2, 0.3, 1.2)
check("anchors_m", anc == {"x_in_m": -0.5, "x_contraction_m": -0.3, "x_throat_m": 0.0, "x_exit_m": 1.2}, str(anc))
s_pts = [-2.5, -2.0, -1.5, -1.0, -0.25, 0.0, 0.5, 1.0, 1.5]
hand = [-0.6, -0.5, -0.4, -0.3, -0.075, 0.0, 0.6, 1.2, 1.8]      # 手計算 (外側は同じ区間の写像の延長)
xm, Tm = WTP.to_x_m({"coord": "s", "coords": s_pts, "Tw_K": [300.0] * len(s_pts)}, anc)
check("s → x [m] が手計算と相対 1e-15 (区切りの点・外側の延長を含む)",
      rel_ok(xm, hand, 1e-15) and xm[5] == 0.0, f"{xm.tolist()}")
anc_d = WTP.anchors_from_design(P_S, D, A)
L_pipe, L_U = float(GEOM["L_pipe"]), float(GEOM["L_U"])
check("anchors_from_design: 区切り = L_pipe・L_U (問題) と壁の終点 x_e (設計) × scale",
      anc_d == {"x_in_m": -(L_pipe * A + L_U * A), "x_contraction_m": -(L_U * A), "x_throat_m": 0.0, "x_exit_m": float(D["wall"].x_e) * A},
      str(anc_d))
xs_d, _ = WTP.to_x_m(ps, anc_d)
check("s の表の区切りの点 (−2・−1・0・1) は区切りの x にそのまま写る",
      xs_d[0] == anc_d["x_in_m"] and xs_d[1] == anc_d["x_contraction_m"] and xs_d[4] == 0.0 and xs_d[-1] == anc_d["x_exit_m"])
check("s の表の写しが試験側の独立な写像と相対 1e-15", rel_ok(xs_d, ref_x_of_table(PROF_S, A, D, GEOM), 1e-15))
ok, msg = raises(lambda: WTP.check_covers([-0.4, 1.3], -0.5, 1.2, "試験の範囲"), "必要な範囲 [-0.5, 1.2] m", "表の範囲 [-0.4, 1.3] m")
check("check_covers: 覆わない表は必要な範囲と表の範囲を示して例外", ok, msg)
ok_cov = True
try:
    WTP.check_covers([-0.5 - 5e-13, 1.2 - 5e-13], -0.5, 1.2, "試験の範囲")
except ValueError:
    ok_cov = False
check("check_covers: 許容 1e-12 m の内側は通す", ok_cov)
check("eval_Tw: 範囲内は np.interp と一致", float(WTP.eval_Tw([0.0, 1.0], [300.0, 400.0], 0.25)) == 325.0)
ok, msg = raises(lambda: WTP.eval_Tw([0.0, 1.0], [300.0, 400.0], [0.5, 1.1]), "表の範囲")
check("eval_Tw: 範囲外は例外 (端値クランプしない)", ok, msg)
ok, msg = raises(lambda: WTP.to_x_m(ps, WTP.anchors_m(0.0, 0.3, 1.2)), "狭義単調増加")
check("to_x_m: 入口直管の長さ 0 では s の表を写せない (例外)", ok, msg)

print("== V1 thermal_bc の Tw_table", flush=True)
for nm, pp_, prof in (("x_m", P_X, PROF_X), ("s", P_S, PROF_S)):
    for tag, sc in (("A", A), ("B", B)):
        tb = np.asarray(WTP.integral_thermal_bc(pp_, D, sc)["Tw_table"])
        x_exp = ref_x_of_table(prof, sc, D, GEOM) / sc
        check(f"Tw_table ({nm}, scale {tag}) の x が x_m / scale と相対 1e-15・T が元の表",
              rel_ok(tb[:, 0], x_exp, 1e-15) and tb[:, 1].tolist() == prof["Tw_K"])
tb_s_A = WTP.integral_thermal_bc(P_S, D, A)["Tw_table"]
tb_s_B = WTP.integral_thermal_bc(P_S, D, B)["Tw_table"]
check("s の表は r_t 単位の Tw_table が scale に依らない (長さのパラスタで伸縮する)", rel_ok(np.asarray(tb_s_A), np.asarray(tb_s_B), 1e-14))

print("== V1 bcond の壁の行", flush=True)
line = P_S.wall_bcond_line(False, phys_id=3, output=1, profile_ok=True)
exp_line = ("{physID: 3, kind: wall_isothermal,  outputHDFflg: 1, ints: {wallProfile: 1}, "
            "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 500.0}}")
check("分布の壁の行: ints {wallProfile: 1}・Ts は代表値", line == exp_line, line)
ld = yaml.safe_load(line)
check("その行が YAML として ints.wallProfile = 1・kind wall_isothermal・Ts 500",
      ld["ints"] == {"wallProfile": 1} and ld["kind"] == "wall_isothermal" and ld["floats"]["Ts"] == 500.0)
ok, msg = raises(lambda: P_S.wall_bcond_line(False, phys_id=3, output=1), "未対応")
check("profile_ok=False (既定; SERN など CSV を書かない runner) は例外", ok, msg)
ok, msg = raises(lambda: _bcond(P_S, euler=False), "未対応")
check("runner_wt._bcond(p, euler=False) の既定も例外", ok, msg)
check("Euler (slip) は分布でも変わらない", P_S.wall_bcond_line(True, phys_id=3, output=1) == P_U.wall_bcond_line(True, phys_id=3, output=1)
      and _bcond(P_S, euler=True) == _bcond(P_U, euler=True))
check("runner_wt._bcond(p, euler=False, profile_ok=True) の壁の行", ("wall:   " + exp_line) in _bcond(P_S, euler=False, profile_ok=True))

print("== V1 CSV と JSON の往復 (write_run_files)", flush=True)
fdir = OUT / "v1_files"
shutil.rmtree(fdir, ignore_errors=True); fdir.mkdir()
files = WTP.write_run_files(fdir, P_X, D, A, wall_id=3)
check("write_run_files の戻り (名前・sha256)", files["json"] == WTP.RUN_JSON and files["csv"] == "wall_profile_3.csv"
      and len(files["json_sha256"]) == 64 and len(files["csv_sha256"]) == 64, str(files)[:120])
txt = (fdir / "wall_profile_3.csv").read_text().splitlines()
hdr = txt[0].split()
rows = np.array([[float(t) for t in ln.split()] for ln in txt[1:]])
rec = json.loads((fdir / WTP.RUN_JSON).read_text())
check("CSV のヘッダは空白区切りの `x Ts` (ソルバの applyBoundaryProfiles は空白で分割する; case/48 run_0015 と同じ)",
      hdr == ["x", "Ts"] and "," not in txt[0])
check("CSV の値を読み直すと JSON の x_m・Tw_x_K と相対 1e-12",
      rel_ok(rows[:, 0], rec["x_m"], 1e-12) and rel_ok(rows[:, 1], rec["Tw_x_K"], 1e-12), f"{rows.shape}")
check("JSON: format・元の表・区切りの点・実寸の表・csv の名前",
      rec["format"] == "forge.wall_temperature_profile/v1" and rec["coord"] == "x_m" and rec["coords"] == PROF_X["x_m"]
      and rec["Tw_K"] == PROF_X["Tw_K"] and rec["anchors_m"] == WTP.anchors_from_design(P_X, D, A)
      and rec["x_m"] == PROF_X["x_m"] and rec["Tw_x_K"] == PROF_X["Tw_K"] and rec["csv"] == "wall_profile_3.csv")
ok, msg = raises(lambda: WTP.write_run_files(fdir, P_U, D, A), "壁温分布")
check("一様の問題で write_run_files は例外", ok, msg)
fdir_s = OUT / "v1_files_s"
shutil.rmtree(fdir_s, ignore_errors=True); fdir_s.mkdir()
WTP.write_run_files(fdir_s, P_S, D, A, wall_id=3)
rows_s = np.array([[float(t) for t in ln.split()] for ln in (fdir_s / "wall_profile_3.csv").read_text().splitlines()[1:]])
check("s の表の CSV の x が試験側の独立な写像と相対 1e-12", rel_ok(rows_s[:, 0], ref_x_of_table(PROF_S, A, D, GEOM), 1e-12)
      and rows_s[:, 1].tolist() == PROF_S["Tw_K"])

# ======================================================================================================
print("== V1b 実効の r_t の A/B (integral_delta_r の scale だけを変える)", flush=True)
init = dict(P_U.raw["deltastar_initializer"])
v1b = {}
for nm, pp_, prof in (("x_m", P_X, PROF_X), ("s", P_S, PROF_S)):
    for tag, sc in (("A", A), ("B", B)):
        res, _drx, info = RA.integral_delta_r(pp_, D, init, scale=sc)
        x_tab = ref_x_of_table(prof, sc, D, GEOM)
        Tref = np.interp(sc * np.asarray(res["x"]), x_tab, np.asarray(prof["Tw_K"], dtype=float))
        md = float(np.max(np.abs(np.asarray(res["Tw"]) - Tref)))
        span = float(np.ptp(res["Tw"]))
        v1b[f"{nm}_{tag}"] = md
        check(f"V1b {nm} {tag} (scale {sc!r} m): max|res[Tw] − 原表(scale·x)| ≤ 1e-8 K", md <= 1e-8,
              f"最大差 {md:.3e} K、温度の幅 {span:.1f} K、{len(res['x'])} 点")
        check(f"V1b {nm} {tag}: 記録の thermal_bc が integral_thermal_bc(scale) と同じ",
              info["thermal_bc"] == WTP.integral_thermal_bc(pp_, D, sc))
        if nm == "x_m" and tag == "B":
            wrong = np.interp(np.asarray(res["x"]), np.asarray(PROF_X["x_m"]) / A, np.asarray(PROF_X["Tw_K"]))
            dw = float(np.max(np.abs(wrong - Tref)))
            v1b["x_m_B_wrong_radius"] = dw
            check("V1b 判別力: 旧半径 (spec.r_throat) で換算する誤実装なら B で 10 K 以上ずれる", dw >= 10.0, f"{dw:.1f} K")
ok, msg = raises(lambda: RA.integral_delta_r(P_XN, D, init, scale=B), "覆わない", "必要な範囲")
check("V1b 負例: x_m の表が B の積分区間を覆わない → 例外 (端値クランプで通らない)", ok, msg)
try:
    resA_n, _, _ = RA.integral_delta_r(P_XN, D, init, scale=A)
    okA = float(np.max(np.abs(np.asarray(resA_n["Tw"]) - np.interp(A * np.asarray(resA_n["x"]), PROF_XN["x_m"], PROF_XN["Tw_K"])))) <= 1e-8
    whyA = ""
except Exception as ex:        # noqa: BLE001
    okA, whyA = False, f"{type(ex).__name__}: {str(ex)[:160]}"
check("V1b 同じ表で A (覆う) は通る", okA, whyA)
sup = DL._sizing_delta_supplier(P_X, init)
_drx, srec = sup(D, B)
tb_sup = np.asarray(srec["settings"]["thermal_bc"]["Tw_table"])
check("寸法の逆算の CFD 前の供給 (_sizing_delta_supplier) も試しの r_t で表を換算する (Tw_table の x = x_m / r_t)",
      rel_ok(tb_sup[:, 0], np.asarray(PROF_X["x_m"]) / B, 1e-15))
SUMMARY["V1b"] = v1b

# ======================================================================================================
print("== V1d 正本の照合 (prepare-only の run)", flush=True)
pf_s = OUT / "c70_prof_s.yaml"
pf_s2 = OUT / "c70_prof_s2.yaml"
runA = OUT / "v1d_run_prof_s"
runB = OUT / "v1d_run_prof_s2"
rcA, logA = prepare_cli(DESIGN, pf_s, runA)
check("V1d: 分布 (s) の問題で --init-integral --prepare-only が通る", rcA == 0, logA[-300:] if rcA else "")
rcB, logB = prepare_cli(DESIGN, pf_s2, runB)
check("V1d: スロートの温度が同じで分布だけ違う表でも通る", rcB == 0, logB[-300:] if rcB else "")
if rcA == 0 and rcB == 0:
    piA = json.loads((runA / "prepare_info.json").read_text())
    wtA = piA["wall_thermal"]
    fA = wtA.get("files", {})
    import hashlib
    sha = lambda f: hashlib.sha256(f.read_bytes()).hexdigest()          # noqa: E731
    check("prepare_info.wall_thermal: mode・代表値・Tw_profile・files (sha256 が実ファイルと一致)",
          wtA.get("mode") == "isothermal" and wtA.get("Tw") == 500.0 and wtA.get("Tw_profile") == ps
          and fA.get("json") == WTP.RUN_JSON and fA.get("csv") == "wall_profile_3.csv"
          and fA.get("json_sha256") == sha(runA / WTP.RUN_JSON) and fA.get("csv_sha256") == sha(runA / "wall_profile_3.csv"),
          str({k: v for k, v in wtA.items() if k != "Tw_profile"})[:200])
    bcA = yaml.safe_load((runA / "bcondConfig.yaml").read_text())
    check("bcondConfig の壁: wall_isothermal + ints {wallProfile: 1}・Ts 500",
          bcA["wall"]["kind"] == "wall_isothermal" and bcA["wall"]["ints"] == {"wallProfile": 1} and bcA["wall"]["floats"]["Ts"] == 500.0)
    drA = json.loads((runA / "delta_r_initial.json").read_text())
    check("積分法 (CFD 前の初期壁) が Tw_table を読んだ記録 (delta_r_initial.json の thermal_bc = integral_thermal_bc(spec.r_throat))",
          drA["thermal_bc"] == WTP.integral_thermal_bc(P_S, D, A))
    recA = WTP.load_run_profile(runA)
    xs_run = np.asarray(recA["x_m"])
    # 最終の物理壁の範囲: 物理壁の x_in〜x_e (prepare の物理壁 = wall_physical.csv は x_throat 以降なので設計壁の端で見る)
    check("load_run_profile: 書き換えていない run は通り、実寸の表が to_x_m と同じ",
          recA is not None and recA["x_m"] == WTP.to_x_m(ps, WTP.anchors_from_design(P_S, D, A))[0].tolist())
    wphys = np.loadtxt(runA / "wall_physical.csv", delimiter=",", skiprows=1)
    check("CSV が物理壁の範囲 (入口端〜壁の終点) を覆う", xs_run[0] <= float(D["wall"].x_in) * A + 1e-12
          and xs_run[-1] >= float(wphys[-1, 0]) - 1e-12, f"表 [{float(xs_run[0])!r}, {float(xs_run[-1])!r}]、物理壁の終点 {float(wphys[-1, 0])!r} m")
    check("一様の既存 run (case/70 run_0010) は None (読むだけ)", WTP.load_run_profile(C70_RUN) is None)

    # --- 負例 1: CSV を 1 値書き換える ---
    t1 = copy_run(runA, "v1d_tamper_csv")
    lines = (t1 / "wall_profile_3.csv").read_text().splitlines()
    xv, tv = lines[3].split()
    lines[3] = f"{xv} {float(tv) + 1.0!r}"
    (t1 / "wall_profile_3.csv").write_text("\n".join(lines) + "\n")
    cfg_before = (t1 / "solverConfig.yaml").read_bytes()
    ok, msg = raises(lambda: WTP.load_run_profile(t1), "CSV が書き換えられた")
    check("負例: CSV を 1 値書き換えた → load_run_profile が止まる", ok, msg)
    rc, log = staged_cli(t1)
    traces = [n for n in LAUNCH_TRACES if (t1 / n).exists()]
    check("負例: 同じ run で run_staged_ns が forge を起動する前に止まる (起動の痕跡なし・solverConfig 不変)",
          rc != 0 and "CSV が書き換えられた" in log and not traces and (t1 / "solverConfig.yaml").read_bytes() == cfg_before,
          f"rc {rc}、痕跡 {traces}、{log.strip().splitlines()[-1][:160] if log.strip() else ''}")

    # --- 負例 2: スロートの温度が同じで分布だけ違う別の run の JSON を写す ---
    t2 = copy_run(runA, "v1d_tamper_json")
    shutil.copy2(runB / WTP.RUN_JSON, t2 / WTP.RUN_JSON)
    recB = json.loads((runB / WTP.RUN_JSON).read_text())
    check("(前提) 2 つの run はスロートの温度が同じで分布が違う",
          recB["Tw_representative_K"] == recA["Tw_representative_K"] == 500.0 and recB["Tw_x_K"] != recA["Tw_x_K"])
    ok, msg = raises(lambda: WTP.load_run_profile(t2), "書き換えられた")
    check("負例: 別の run の JSON を写した → load_run_profile が止まる", ok, msg)
    rc, log = staged_cli(t2)
    traces = [n for n in LAUNCH_TRACES if (t2 / n).exists()]
    check("負例: その run で run_staged_ns が forge を起動する前に止まる", rc != 0 and "書き換えられた" in log and not traces,
          f"rc {rc}、痕跡 {traces}")
    t2b = copy_run(runA, "v1d_tamper_json_csv")
    shutil.copy2(runB / WTP.RUN_JSON, t2b / WTP.RUN_JSON)
    shutil.copy2(runB / "wall_profile_3.csv", t2b / "wall_profile_3.csv")
    ok, msg = raises(lambda: WTP.load_run_profile(t2b), "prepare_info と合わない")
    check("負例: 別の run の JSON と CSV を両方写した → 止まる (prepare_info の sha256 と合わない)", ok, msg)

    # --- 負例 3: 欠け・bcond の不一致・記録の無い分布 ---
    t3 = copy_run(runA, "v1d_missing_csv")
    (t3 / "wall_profile_3.csv").unlink()
    ok, msg = raises(lambda: WTP.load_run_profile(t3), "欠け")
    check("負例: CSV が無い → 止まる", ok, msg)
    t4 = copy_run(runA, "v1d_bcond_noflag")
    (t4 / "bcondConfig.yaml").write_text((t4 / "bcondConfig.yaml").read_text().replace("ints: {wallProfile: 1}", "ints: "))
    ok, msg = raises(lambda: WTP.load_run_profile(t4), "wallProfile")
    check("負例: bcond から wallProfile を消した → 止まる (ソルバが CSV を読まない)", ok, msg)
    t5 = OUT / "v1d_stray_json"
    shutil.rmtree(t5, ignore_errors=True); t5.mkdir()
    (t5 / "prepare_info.json").write_text(json.dumps({"wall_thermal": {"mode": "adiabatic"}}))
    shutil.copy2(runA / WTP.RUN_JSON, t5 / WTP.RUN_JSON)
    ok, msg = raises(lambda: WTP.load_run_profile(t5), "記録が無い")
    check("負例: 記録 (files) の無い run に壁温の JSON がある → 止まる", ok, msg)

    # --- 段の間 (restart の後) の書き換え: 各段の起動前に照合する (forge と段の検査は偽物に差し替える) ---
    t6 = copy_run(runA, "v1d_tamper_between_stages")
    calls = []

    def fake_run_forge(run_dir):
        calls.append(Path(run_dir).name)
        (Path(run_dir) / "res_3000.h5").write_bytes(b"")
        if len(calls) == 1:          # 1 段目の後 (次の段の restart の前後) に CSV を書き換える
            cp = Path(run_dir) / "wall_profile_3.csv"
            cp.write_text(cp.read_text().replace("500.0", "501.0"))
        return 0
    saved = (RA.run_forge, RA.stage_gate, RA._restart_same_mesh)
    RA.run_forge, RA.stage_gate, RA._restart_same_mesh = fake_run_forge, (lambda *a, **k: []), (lambda *a, **k: None)
    try:
        ok, msg = raises(lambda: RA.run_staged_ns(t6, stages="full"), "書き換えられた")
    finally:
        RA.run_forge, RA.stage_gate, RA._restart_same_mesh = saved
    check("段の間に CSV を書き換えた → 次の段の起動前に止まる (forge の呼び出しは 1 回だけ)", ok and len(calls) == 1,
          f"呼び出し {len(calls)} 回; {msg}")

    # --- 対照: 書き換えていない run は照合を通り、起動 (run_case.sh) まで進む ---
    t7 = copy_run(runA, "v1d_control_launch")
    rc, log = staged_cli(t7, stages="none")
    check("対照: 書き換えていない run は照合を通り起動まで進む (FORGE_BIN が無いので forge は失敗し、起動の痕跡が残る)",
          (t7 / "run_case_stdout.log").exists() and "書き換えられた" not in log and "RUN_STAGED_RC" in log,
          log.strip().splitlines()[-1][:160] if log.strip() else "")

# ======================================================================================================
print("== V4 後方互換 (一様・断熱は変更前と同一)", flush=True)


def ref_source(rel: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), "show", f"{REF_COMMIT}:{rel}"], capture_output=True, text=True, check=True).stdout


ref_pd = types.ModuleType("forge_design._ref_probdef")
ref_pd.__package__ = "forge_design"
sys.modules[ref_pd.__name__] = ref_pd
exec(compile(ref_source("design/forge_design/probdef.py"), f"{REF_COMMIT[:8]}:probdef.py", "exec"), ref_pd.__dict__)
_tree = ast.parse(ref_source("design/forge_design/evaluate/runner_wt.py"))
_fn = next(n for n in _tree.body if isinstance(n, ast.FunctionDef) and n.name == "_bcond")
_ns = {"Problem": object}
exec(compile(ast.Module(body=[_fn], type_ignores=[]), f"{REF_COMMIT[:8]}:runner_wt._bcond", "exec"), _ns)
ref_bcond = _ns["_bcond"]
pf45 = write_problem(C45, yaml.safe_load(C45.read_text())["spec"].get("wall_thermal") or {"mode": "adiabatic"}, "c45_prod_abs", abs_pin=True)
raw45 = yaml.safe_load(C45.read_text())
if "wall_thermal" not in raw45["spec"]:
    # 元の YAML に wall_thermal が無い (既定 = 断熱) なら写しからも消す (キーの有無まで同じにする)
    r_ = yaml.safe_load(pf45.read_text()); r_["spec"].pop("wall_thermal", None)
    pf45.write_text(yaml.safe_dump(r_, sort_keys=False, allow_unicode=True))
v4 = {}
for label, path in (("case/70 c2b (一様 350 K)", C70), ("case/45 生産 (断熱)", C45)):
    pn, po = load_problem(path), ref_pd.load_problem(path)
    same_wt = pn.wall_thermal == po.wall_thermal
    lines_ok = all(pn.wall_bcond_line(eu, phys_id=3, output=1) == po.wall_bcond_line(eu, phys_id=3, output=1)
                   == pn.wall_bcond_line(eu, phys_id=3, output=1, profile_ok=True) for eu in (False, True))
    bc_ok = (_bcond(pn, euler=False) == ref_bcond(po, False) == _bcond(pn, euler=False, profile_ok=True)
             and _bcond(pn, euler=True) == ref_bcond(po, True))
    check(f"V4 {label}: wall_thermal・壁の行・_bcond(p, euler=False) の文字列が変更前と同一", same_wt and lines_ok and bc_ok,
          f"{pn.wall_thermal}")
    if label.startswith("case/45"):
        pd_ = load_problem(pf45)
        d_ = RA.design_chain(pd_) if (C45_PIN_DIR / raw45["geometry"]["initial_line_run"]).exists() else None
    else:
        d_ = D
    sc0 = float(pn.spec["r_throat"])
    tbc_ok = True
    for sc in (sc0, 1.2 * sc0):
        new = WTP.integral_thermal_bc(pn, d_, sc)
        old = po.wall_thermal_bc_integral
        tbc_ok &= (new == old and json.dumps(new) == json.dumps(old) and type(new) is dict)
    check(f"V4 {label}: integral_thermal_bc の dict が変更前の wall_thermal_bc_integral と同一 (scale 1 倍・1.2 倍)",
          tbc_ok, f"{po.wall_thermal_bc_integral}{'' if d_ is not None else ' (design_chain は CFD ピンの run が無いので省略)'}")
    v4[label] = bool(same_wt and lines_ok and bc_ok and tbc_ok)
# 既存の一様の run の bcond とバイト同一 (case/70 run_0010 は同じ問題 YAML から作った run)
if (C70_RUN / "bcondConfig.yaml").exists():
    bc_now = RA._bcond_with_species(_bcond(P_U, euler=False, profile_ok=True), RA._tp_species_Y(P_U))
    check("V4 case/70: 今の prepare_ns の bcond 文字列が既存の run_0010 の bcondConfig.yaml とバイト同一",
          bc_now.encode() == (C70_RUN / "bcondConfig.yaml").read_bytes())

# 変更前のコードと今のコードで同じ問題を prepare-only し、生成物がバイト同一
ref_root = OUT / "ref_code"
shutil.rmtree(ref_root, ignore_errors=True); ref_root.mkdir()
# 設計側の種データ (gas/semiperfect.py が design の 1 つ上の solver_density_cuda/data を読む) も同じ commit から
arc = subprocess.run(f"git -C {REPO} archive {REF_COMMIT} design solver_density_cuda/data | tar -x -C {ref_root}",
                     shell=True, capture_output=True, text=True)
check("V4: 変更前のコード (git archive REF_COMMIT design solver_density_cuda/data) を展開", arc.returncode == 0, arc.stderr[-200:])
CMP_FILES = ("bcondConfig.yaml", "solverConfig.yaml", "delta_r_initial.json", "delta_r_initial.csv", "wall_design.csv",
             "wall_physical.csv", "target_axis_M.csv", "wall_repr.json", "probe.yaml", "MESH_QUALITY.txt")
def _cmp_bytes(f, path):
    # NS の本段の step 数 (2026-10-14 に既定を 60,000 に変えた、runner_axismach.NS_STEPS_DEFAULT) は比べない — 下で値を別に確かめる
    b = path.read_bytes()
    return re.sub(rb"nStepOuter: \d+", b"nStepOuter: N", b) if f == "solverConfig.yaml" else b


for label, prob in (("case/70 c2b (一様 350 K)", C70), ("case/45 生産 (断熱)", pf45)):
    if label.startswith("case/45") and not (C45_PIN_DIR / raw45["geometry"]["initial_line_run"]).exists():
        print(f"   {label}: CFD ピンの run が無いので prepare の比較を飛ばした")
        continue
    tag = "c70" if label.startswith("case/70") else "c45"
    r_old, r_new = OUT / f"v4_{tag}_before", OUT / f"v4_{tag}_after"
    rc_o, log_o = prepare_cli(ref_root / "design", prob, r_old)
    # 今のコードは出口の切り方を直す前の規則に戻して比べる (§5.1 #4c: c2b・case/45 生産は直す前の切り方が x_F の壁点を落としていた;
    # case/70 exit_cut_ab.py)。今の切り方の差は下で別に見る
    rc_n, log_n = prepare_cli(DESIGN, prob, r_new, old_exit_cut=True)
    if rc_o or rc_n:
        check(f"V4 {label}: prepare-only (変更前・今) が通る", False, (log_o if rc_o else log_n)[-300:])
        continue
    diff = [f for f in CMP_FILES if (r_old / f).exists() != (r_new / f).exists()
            or ((r_old / f).exists() and _cmp_bytes(f, r_old / f) != _cmp_bytes(f, r_new / f))]
    pi_o = json.loads((r_old / "prepare_info.json").read_text()); pi_n = json.loads((r_new / "prepare_info.json").read_text())
    for pi_ in (pi_o, pi_n):          # 実行環境 (git の状態・時刻) と壁の当てはめの所要時間 [s] (壁時計) は比べない
        pi_.pop("environment", None)
        if isinstance(pi_.get("physical_wall"), dict) and isinstance(pi_["physical_wall"].get("fit"), dict):
            pi_["physical_wall"]["fit"].pop("build_seconds", None)
    # mesh 欄の節点数と上限 (2026-10-14 に足したキー、plan tooling-nozzle-mesh-size-guard §4.2) は変更前に無いので比べない。
    # 節点数は今の記録が実際の座標数であることだけを見る
    nodes_ok = pi_n.get("mesh", {}).get("n_nodes") == pi_n.get("mesh", {}).get("ni", 0) * pi_n.get("mesh", {}).get("nj", 0)
    for k in ("n_nodes", "max_nodes", "max_nodes_reason", "max_nodes_source"):
        pi_n.get("mesh", {}).pop(k, None)
    # NS の本段の step 数: 変更前は evaluate.nStepOuter (雛形の 12,000) を読み、今は既定 60,000 (evaluate.ns_steps が無いとき)
    ns_o, ns_n, ns_src = pi_o.pop("nStepOuter", None), pi_n.pop("nStepOuter", None), pi_n.pop("nStepOuter_source", None)
    check(f"V4 {label}: NS の本段の step 数は既定 {RA.NS_STEPS_DEFAULT} (evaluate.ns_steps が無いとき; 変更前は evaluate.nStepOuter)",
          ns_n == RA.NS_STEPS_DEFAULT and str(ns_src).startswith("既定"), f"変更前 {ns_o} → 今 {ns_n} ({ns_src})")
    # 停止規則の設定 (2026-10-14 に足したキー、plan tooling-nozzle-cooled-wall-ns-settling §4.6) は変更前に無いので比べない。
    # 既定 (要求あり・上限 200,000・使わないときの本段 = NS_STEPS_DEFAULT) であることを別に見る
    sr_n = pi_n.pop("stop_rule", None) or {}
    check(f"V4 {label}: 停止規則の設定は既定 (要求あり・上限 200,000・使わないときの本段 {RA.NS_STEPS_DEFAULT})",
          sr_n.get("requested") is True and sr_n.get("cap") == 200000 and sr_n.get("fallback_steps") == RA.NS_STEPS_DEFAULT
          and "stop_rule" not in pi_o, str({k: sr_n.get(k) for k in ("requested", "use", "cap", "cap_source", "fallback_steps", "reason")}))
    pi_same = nodes_ok and json.dumps(pi_o, sort_keys=True) == json.dumps(pi_n, sort_keys=True)
    pi_diff = sorted(k for k in set(pi_o) | set(pi_n) if pi_o.get(k) != pi_n.get(k))
    extra = [f.name for f in r_new.iterdir() if f.name == WTP.RUN_JSON or f.name.startswith("wall_profile_")]
    present = [f for f in CMP_FILES if (r_new / f).exists()]
    check(f"V4 {label}: prepare-only の生成物が変更前とバイト同一 ({len(present)} 個 + prepare_info から environment・"
          "physical_wall.fit.build_seconds・mesh の節点数と上限・stop_rule を除いたもの)、壁温の JSON・CSV を書かない", not diff and pi_same and not extra,
          f"違うファイル {diff}、prepare_info の違うキー {pi_diff}、余分 {extra}、mesh.n_nodes = ni × nj {nodes_ok}")
    check(f"V4 {label}: 一様・断熱の run の load_run_profile は None", WTP.load_run_profile(r_new) is None)
    v4[label + " prepare"] = bool(not diff and pi_same and not extra)
    # 今の切り方 (§5.1 #4c): 違うのは壁の長さに依るファイルだけで、設計壁の終点が直す前より下流 (x_F) に動く
    r_cut = OUT / f"v4_{tag}_after_exitcut"
    rc_c, log_c = prepare_cli(DESIGN, prob, r_cut)
    if rc_c:
        check(f"V4 {label}: 今の切り方の prepare-only が通る", False, log_c[-300:])
        continue
    wall_dep = {"delta_r_initial.json", "delta_r_initial.csv", "wall_design.csv", "wall_physical.csv", "wall_repr.json", "MESH_QUALITY.txt"}
    diff_c = [f for f in CMP_FILES if (r_old / f).exists() != (r_cut / f).exists()
              or ((r_old / f).exists() and _cmp_bytes(f, r_old / f) != _cmp_bytes(f, r_cut / f))]
    dom_o = json.loads((r_old / "wall_repr.json").read_text())["domain"]
    dom_c = json.loads((r_cut / "wall_repr.json").read_text())["domain"]
    check(f"V4 {label}: 今の切り方は壁の長さに依るファイルだけが違い、設計壁の終点が直す前より下流 (入口端は同じ)",
          set(diff_c) <= wall_dep and dom_c[0] == dom_o[0] and dom_c[1] > dom_o[1],
          f"違うファイル {diff_c}、終点 {dom_o[1]!r} → {dom_c[1]!r} (差 {dom_c[1] - dom_o[1]:.6e} r_t)")
SUMMARY["V4"] = v4

# --- codex result レビュー 2026-10-13 m1: s = 0 を点として含まない s の表は拒否 (代表値と実際のスロートの温度の食い違いを作らない) ---
r_, msg_ = raises(lambda: WTP.parse_profile({"s": [-2.0, -1.0, 1.0], "Tw_K": [300.0, 400.0, 600.0]}), "s = 0 を点として含む")
check("m1: s = 0 を含まない s の表 (codex の再現例 s [-2,-1,1]・T [300,400,600]) は拒否", r_, msg_)
p0_ = WTP.parse_profile({"s": [-2.0, -1.0, 0.0, 1.0], "Tw_K": [300.0, 400.0, 440.0, 600.0]})
anc_ = WTP.anchors_m(0.2, 0.3, 1.2)
xm_, tm_ = WTP.to_x_m(p0_, anc_)
check("m1: s = 0 を含む表では代表値 = 実寸に写した表の x = 0 の値", WTP.representative_Tw(p0_) == float(np.interp(0.0, xm_, tm_)),
      f"代表値 {WTP.representative_Tw(p0_)!r}・x = 0 {float(np.interp(0.0, xm_, tm_))!r}")

print("\n要約:", json.dumps(SUMMARY, ensure_ascii=False, indent=1))
print(f"\n{'ALL OK' if FAIL == 0 else f'{FAIL} FAIL'}")
if FAIL == 0 and not KEEP_OUT:
    shutil.rmtree(OUT, ignore_errors=True)        # 一時ディレクトリは合格なら消す (不合格なら診断用に残す)
else:
    print(f"出力を残した: {OUT}")
sys.exit(1 if FAIL else 0)
