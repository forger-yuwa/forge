"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36b・#36c (§4.9 codex diagnose 2)・#36e (codex diagnose 3、§6 V-h370): 370 K・H・B0 の NS の
問題と、H の固定の Euler の参照 (補正 δ_r = 0) の問題を作り、読み込んだ辞書の差が登録の差だけであることを検査し、CFD 0 step (変換器・forge を
使わない) で準備の中身を確かめる。

問題 (make — 既にあれば生成した中身と一致することを確かめ、違えば止める; 手で編集しない):
- NS: problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml ← problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml (#34、H・Gref・300 K) の行単位の複製。
  差は name・spec.wall_thermal.Tw (300 → 370)・evaluate (cfl_main を外し、implicit_relax 0.7・main_line (本段だけ B0)・nStepOuter 200000・
  outStepInterval 2500・停止規則のキー (#36e: stop_rule・stop_rule_cap 200000・stop_rule_hold_steps 40000・stop_rule_thresholds; §6 V-h370 の値))。
  evaluate は停止規則のキーを除いて V-g0c-chain (problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml) と V-tw370 の C50 の evaluate と同じ
  (リミッタの基準値は今の形 C の系譜の値に固定; codex diagnose 3 で決着)。
- Euler: problem_d155_euler_h_ref.yaml ← run_0174_euler_v5d_M_r1 を作った問題 problem_d155_euler_v5d_M.yaml。この作業ツリーには無いので
  moc_v5d.problem_texts (make-problems と同じ関数、--md-offset 6.8825162455159465e-06) で作り直し、V5d の起動の記録
  (_band_ab/moc_v5d_run.log の input 行) の sha256 cfee7f05… と一致することを確かめてから使う。差は name・spec.r_throat (NS と同じ値)・
  geometry.physical_upstream (NS と同じ H)・mesh_euler.throat_refine (4 → 5; #36e、ni 2000 × nj 97 は不変)・evaluate.nStepOuter (12000 → 60000;
  #36e、本段の上限。出力 1000・cfl_main 2 は元のまま)。r_throat を NS に揃えるのは、抽出器の流量の比 (massflow_ratio) が実寸で取るため
  (r_t 単位の Euler の解は r_t によらない)。δ_r = 0 は YAML のキーではなく Euler の prepare の経路の定義 (runner_axismach.euler_reference_wall、
  #36e の案 A: H の上流を付けた設計壁)。
δ_r (NS): prepare_ns(delta_r_csv=_band_ab/half_inlet/delta_r_initial_C.csv) の引数で渡す (YAML のキーにしない; prepare_ns の既存の契約、
  prepare_info.json の delta_r_csv に残る)。sha256 24251cb6… を検査する。

check: 上の差の検査 (読み込みは重複キーを拒む yaml_strict)、NS の evaluate が V-g0c-chain・V-tw370 C50 と同じ、Euler と NS の設計の入力
  (gas・spec (wall_thermal 以外)・dv・geometry (NS だけの pw_upstream・physical_wall_repr 以外)) が同じ、main_line の解決 (`_main_line_spec`)。
prep-check (CFD 0 step): prepare_ns の変換の前までを同じ関数で再現する (design_chain・CFD ピンの照合・build_physical_wall・mesh_params・
  generate_axisym_mesh・_bcond・_ns_configs・ns_stage_configs)。格子の品質 (check_mesh_quality.metrics_vectorized を生成の座標に、
  ar_max 5000)・第一層の距離・S1・S2・本段の実効設定・bcond の壁温・300 K の H (#34) と格子がビット同一・生成器の sha256
  (prepare_info.mesh.hashes.generator_sha256 と同じ式) を記録する。Euler は design_chain の設計壁 S が NS と同じか、Euler の参照の壁
  (euler_reference_wall、#36e の案 A の正式の経路) のゲート、診断の参照の壁 (#34 の A の定義 HRefWall) との一致、壁ファイルの保存 → 復元で
  格子の座標がビット一致、Euler の格子 (mesh_euler、throat_refine 5) の品質 (生成の座標)、抽出器の δ_in と NS の壁の入力 δ_r の照合
  (metrics.delta_in_check、x ≥ 0 で ≤ 0.001 mm・全列の被覆)。NS の格子で一定 x の断面への補間の端の補間範囲の不足の列の数 (受け渡しの抽出の
  remap_edge="missing" で欠測になる列; 0 であること) と停止規則の要求の断面の位置 (x_fracs × x_w) も記録する。
  変換後の品質 (check_mesh_quality --ar-max 5000)・変換器の wall_dist (局所の法線の距離)・初期場は AWS で prepare_ns・prepare を回して確かめる。

usage: /home/sano/work/forge/design/.venv-opt/bin/python h370_problems.py make|check|prep-check [--pin <run_0062 のディレクトリ>]
  --pin: CFD ピンの凍結元 (geometry.initial_line_run)。無ければこの case の run_0062_euler_wallfit_fit_r1_ext6k。照合は消さない
  (読み込む写しの initial_line_run だけを絶対パスに書き換える)。
出力: _band_ab/half_inlet/h370_problems.json (check)・h370_prep_check_tr5.json と h370_delta_in_tr5.csv (prep-check; #36e。throat_refine 4 の
  旧い FAIL の記録 h370_prep_check.json・h370_delta_in.csv は消さずに残す)。合格で rc 0、不成立で rc 1。"""
import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))

H300 = HERE / "problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml"              # #34 の H (300 K・cfl_main 1・B0 なし)
NS = HERE / "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"             # #36c
EU = HERE / "problem_d155_euler_h_ref.yaml"                                   # #36b
G0C_CHAIN = HERE / "problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml"       # V-g0c-chain (B0 の evaluate の登録の元)
C50 = HERE / "problem_d155_ns_prod_coldmesh_tw370_g0c_c50.yaml"               # V-tw370 の C50
CSV = HERE / "_band_ab" / "half_inlet" / "delta_r_initial_C.csv"
CSV_SHA = "24251cb647226d6523109655c08c05210fdb8ba4a669aeaecb4407202945b918"
V5D_MD = "6.8825162455159465e-06"
V5D_M_SHA = "cfee7f057648c21ef1cfd51a167397c4490ab2ee6a6f8fee31df7fcf84fdaf92"   # _band_ab/moc_v5d_run.log の input 行 (run_0174 の問題)
V5D_LOG = HERE / "_band_ab" / "moc_v5d_run.log"
OUT = HERE / "_band_ab" / "half_inlet"
PIN_DEFAULT = HERE / "run_0062_euler_wallfit_fit_r1_ext6k"
TW_NEW = 370.0
MISSING = "<無し>"
# 停止規則 (§6 V-h370 の登録の値。キー名は stop_rule.py の thresholds・check_evaluate_keys。既定と同じ値も明示する; §5.1 #36e (3))
STOP_RULE_X_FRACS = [0.419738616669543, 0.734542579171700, 0.986385749173426]     # 40・70・94 r_t / x_w (95.2974)
STOP_RULE_THRESHOLDS = {"window": 20000, "min_points": 9, "r_mean_max": 5.0e-4, "r_amp_max": 5.0e-3, "theta_drift_max": 5.0e-4,
                        "theta_amp_max": 5.0e-3, "x_fracs": STOP_RULE_X_FRACS}
STOP_RULE_LINES = ["  # 本段の停止規則 (§6 V-h370 の登録の値; §5.1 #36e)。キー名は stop_rule.py の thresholds・check_evaluate_keys、既定と同じ値も明示",
                   "  stop_rule: true",
                   "  stop_rule_cap: 200000   # 本段の上限 (nStepOuter と同じ。自動の延長はしない)",
                   "  stop_rule_hold_steps: 40000   # 最初の到達の後 4 万 step すべての出力で保つ (20 万の内数、外れたら数え直す)",
                   "  stop_rule_thresholds: {window: 20000, min_points: 9, r_mean_max: 5.0e-4, r_amp_max: 5.0e-3, theta_drift_max: 5.0e-4, "
                   "theta_amp_max: 5.0e-3, x_fracs: [0.419738616669543, 0.734542579171700, 0.986385749173426]}   # x_fracs = 40・70・94 r_t / x_w 95.2974"]
STOP_RULE_KEYS = ("stop_rule", "stop_rule_cap", "stop_rule_hold_steps", "stop_rule_thresholds")
MAIN_LINE_LINE = ("  main_line: {line: dir, cap: 50, itj: 5, cfl: 4, relax: 0.7, limiter_ref: {ro: 0.8739869154, p: 387491.0573, a: 359.7712185}, "
                  "extra_fields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}")
EVAL_NEW = ["evaluate:",
            "  nStepOuter: 200000   # 本段の step 数の上限 (S1・S2 の 3,000 step ずつは別枠、§6 V-h370 の予算。停止規則で止める)",
            "  outStepInterval: 2500   # 停止規則の出力の間隔",
            "  implicit_relax: 0.7   # main_line.relax と同じ値 (前の段 S1・S2 にも効く; 生産の系譜 run_0182 の cold_pair.py prep と同じ)",
            "  # 本段だけ B0 (cold_cfl.py prep --line dir --itj 5 --cap 50 --cfl 4 --limiter-ref-from … --extra … と同じキー・書き方; §4.8-2・§5.1 #30)。",
            "  # リミッタの基準値は今の形 C の系譜 (V-g0c-chain・V-tw370 の C50) の値に固定 (codex diagnose 3 で決着: 今回の全本段で C 由来の値; §4.9)",
            MAIN_LINE_LINE] + STOP_RULE_LINES
NS_ALLOWED = {"spec.wall_thermal.Tw", "evaluate.cfl_main", "evaluate.implicit_relax", "evaluate.main_line", "evaluate.nStepOuter",
              "evaluate.outStepInterval"} | {f"evaluate.{k}" for k in STOP_RULE_KEYS}
EU_ALLOWED = {"spec.r_throat", "geometry.physical_upstream", "mesh_euler.throat_refine", "evaluate.nStepOuter"}
EU_THROAT_REFINE = 5.0          # §5.1 #36e (codex diagnose 3): 4 → 5 で生成の座標の max|δ_in − δ_r| 1.137e-3 → 7.96e-4 mm (codex の再計算)
EU_NSTEP = 60000                # §6 V-h370: Euler の本段の上限 (出力 1000・cfl_main 2 は元の問題のまま)
NS_ONLY_GEOM = {"pw_upstream", "physical_wall_repr"}


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sha_text(t: str) -> str:
    return hashlib.sha256(t.encode()).hexdigest()


def load_strict(text: str) -> dict:
    import yaml_strict
    return yaml_strict.load(text)


def diff_paths(a, b, path=()):
    """2 つの辞書の差 {パス: (a の値, b の値)}。片側にしか無いキーは MISSING (型が違うのも差; tw370_yaml_check.py と同じ定義)。"""
    out = {}
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            out.update(diff_paths(a.get(k, MISSING), b.get(k, MISSING), path + (k,)))
        return out
    if type(a) is not type(b) or a != b:
        out[path] = (a, b)
    return out


def replace_line(lines, prefix, new):
    idx = [i for i, ln in enumerate(lines) if ln.startswith(prefix)]
    if len(idx) != 1:
        raise SystemExit(f"{prefix!r} で始まる行が {len(idx)} 行 — 止める")
    lines[idx[0]] = new
    return idx[0]


def v5d_m_text() -> str:
    """run_0174 の問題 (problem_d155_euler_v5d_M.yaml) の本文を moc_v5d.problem_texts で作り直し、V5d の起動の記録の sha256 と照合する。"""
    sys.path.insert(0, str(HERE))
    import moc_v5d
    t = moc_v5d.problem_texts(V5D_MD)["M"]
    if sha_text(t) != V5D_M_SHA:
        raise SystemExit(f"作り直した problem_d155_euler_v5d_M.yaml の sha256 {sha_text(t)[:16]}… が V5d の記録 {V5D_M_SHA[:16]}… と違う — 止める")
    if f"problem_d155_euler_v5d_M.yaml {V5D_M_SHA}" not in V5D_LOG.read_text():
        raise SystemExit(f"{V5D_LOG.name} に problem_d155_euler_v5d_M.yaml {V5D_M_SHA[:16]}… の行が無い — 止める")
    return t


def ns_text() -> str:
    src = H300.read_text()
    lines = src.splitlines()
    replace_line(lines, "name: ", "name: problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0")
    replace_line(lines, "  wall_thermal: ", f"  wall_thermal: {{mode: isothermal, Tw: {TW_NEW}}}   # 370 K (ユーザ決定 2026-10-13、V-tw370 と同じ壁温)")
    i0 = replace_line(lines, "evaluate:", "evaluate:")
    old = lines[i0 + 1:i0 + 4]
    if old != ["  nStepOuter: 12000", "  outStepInterval: 1000", "  cfl_main: 1.0"]:
        raise SystemExit(f"{H300.name} の evaluate の先頭 3 行が想定と違う ({old}) — 止める")
    lines[i0:i0 + 4] = EVAL_NEW
    hdr = ("# plan tooling-nozzle-isothermal-wall-chain §4.9 (codex diagnose 2・3)・§5.1 #36c・#36e: 370 K・H・Gref・B0 の最初の固定壁の pass の問題\n"
           "#   (H・B0 による暫定の設計、§6 V-h370)。h370_problems.py make が problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml (#34) から作る (手で\n"
           "#   編集しない)。元との差は name・spec.wall_thermal.Tw (300 → 370)・evaluate (cfl_main を外し、implicit_relax 0.7・main_line (本段だけ B0)・\n"
           "#   nStepOuter 200000・outStepInterval 2500・停止規則のキー (#36e)) だけ。evaluate は停止規則のキーを除いて V-g0c-chain・V-tw370 の C50 の\n"
           "#   evaluate と同じ。h370_problems.py check が読み込んだ辞書の差で確かめる。\n"
           f"#   δ_r は prepare_ns(delta_r_csv=_band_ab/half_inlet/delta_r_initial_C.csv) の引数で渡す (sha256 {CSV_SHA[:16]}…、#34 で作り直した断熱の\n"
           "#   積分法の表。「断熱の積分法から作った固定の形に 370 K を課す試験」; deltastar_initializer は delta_r_csv があると使われない)。\n"
           "#   壁温は bcondConfig の壁の行 (wall_isothermal Ts) にだけ出る。\n"
           "# ---- 以下は元の YAML ----\n")
    return hdr + "\n".join(lines) + "\n"


def eu_text(ns_doc: dict) -> str:
    src = v5d_m_text()
    lines = src.splitlines()
    # mesh_euler.throat_refine 4 → 5 (#36e)。上の mesh (NS 用 G1) にも同じ前置きの行があるので mesh_euler の区画の中の行だけを置き換える
    i_me = [i for i, ln in enumerate(lines) if ln.startswith("mesh_euler:")]
    if len(i_me) != 1:
        raise SystemExit(f"mesh_euler: の行が {len(i_me)} 行 — 止める")
    j = i_me[0] + 1
    while j < len(lines) and lines[j].startswith("  "):
        if lines[j].startswith("  throat_refine: "):
            if lines[j] != "  throat_refine: 4.0":
                raise SystemExit(f"mesh_euler.throat_refine の行が想定と違う ({lines[j]!r}) — 止める")
            lines[j] = (f"  throat_refine: {EU_THROAT_REFINE}   # 4 → 5 (§5.1 #36e、codex diagnose 3: スロートの Euler の列を密にして生成の座標の "
                        "max|δ_in − δ_r| 1.137e-3 → 7.96e-4 mm ≤ 0.001 mm。ni・nj は不変)")
            break
        j += 1
    else:
        raise SystemExit("mesh_euler の区画に throat_refine が無い — 止める")
    i_ev = replace_line(lines, "evaluate:", "evaluate:")
    if lines[i_ev + 1:i_ev + 4] != ["  nStepOuter: 12000", "  outStepInterval: 1000", "  cfl_main: 2.0"]:
        raise SystemExit(f"Euler の evaluate の先頭 3 行が想定と違う ({lines[i_ev + 1:i_ev + 4]}) — 止める")
    lines[i_ev + 1] = f"  nStepOuter: {EU_NSTEP}   # 本段の上限 (§6 V-h370、#36e)。出力 1000・cfl_main 2 は元のまま"
    replace_line(lines, "name: ", "name: isobutane_m6_d155_euler_h_ref")   # 既存の Euler の問題の流儀 (isobutane_m6_d155_euler_*)
    rt = ns_doc["spec"]["r_throat"]
    replace_line(lines, "  r_throat: ", f"  r_throat: {rt!r}      # [m] NS (problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml) と同じ値 (抽出器の流量の比は実寸で取る)")
    pu = ns_doc["geometry"]["physical_upstream"]
    i = replace_line(lines, "  wall_fit_mono_r2: ", [ln for ln in lines if ln.startswith("  wall_fit_mono_r2: ")][0])
    lines.insert(i + 1, f"  physical_upstream: {{r_inlet: {pu['r_inlet']}, L_U: {pu['L_U']}, L_pipe: {pu['L_pipe']}}}   # H (NS と同じ)。"
                        "Euler の prepare は補正 δ_r ≡ 0 の物理壁で格子を作る (euler_reference_wall)。r_inlet・L_U・L_pipe (上) は設計壁と CFD ピンの照合の値")
    hdr = ("# plan tooling-nozzle-isothermal-wall-chain §4.9 (codex diagnose 2・3)・§5.1 #36b・#36e: H の固定の Euler の参照 (補正 δ_r = 0、§6 V-h370)。\n"
           "#   h370_problems.py make が run_0174_euler_v5d_M_r1 の問題 problem_d155_euler_v5d_M.yaml (moc_v5d.py make-problems --md-offset\n"
           f"#   {V5D_MD} の本文、sha256 {V5D_M_SHA[:16]}… を確かめて使う) から作る (手で編集しない)。元との差は name・spec.r_throat (NS と同じ)・\n"
           "#   geometry.physical_upstream (NS と同じ H)・mesh_euler.throat_refine (4 → 5)・evaluate.nStepOuter (12000 → 60000) だけ。δ_r = 0 は YAML の\n"
           "#   キーではなく Euler の prepare の経路の定義 (案 A、H の上流を付けた設計壁): x ≥ 0 は凍結した設計壁 S (main の出口の終点の修正の後)、\n"
           "#   x < 0 は直管 r 3.2425 と (−9, 3.2425, 0, 0) → x = 0 の S の (r, r′, r″) の 5 次 Hermite。格子は mesh_euler。\n"
           "#   CFD ピンの照合 (initial_line: cfd・run_0062) は geometry の r_inlet・L_U・L_pipe (凍結した設計) のまま。\n"
           "# ---- 以下は元の YAML ----\n")
    return hdr + "\n".join(lines) + "\n"


def write_or_check(path: Path, text: str) -> str:
    if path.exists():
        if path.read_text() != text:
            raise SystemExit(f"{path.name} が既にあり、生成した中身と違う — 止める (消してから回すか、差を確かめる)")
        return "unchanged"
    path.write_text(text)
    return "written"


def cmd_make() -> int:
    nt = ns_text()
    ns_doc = load_strict(nt)
    et = eu_text(ns_doc)
    print(json.dumps({"ns": {"file": NS.name, "status": write_or_check(NS, nt)}, "euler": {"file": EU.name, "status": write_or_check(EU, et)}},
                     ensure_ascii=False))
    return cmd_check()


def cmd_check(out_json=None) -> int:
    fails, rec = [], {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36b・#36c・#36e", "tool": "h370_problems.py check"}
    D = {"H300": load_strict(H300.read_text()), "NS": load_strict(NS.read_text()), "EU": load_strict(EU.read_text()),
         "V5D_M": load_strict(v5d_m_text()), "G0C_CHAIN": load_strict(G0C_CHAIN.read_text()), "C50": load_strict(C50.read_text())}
    rec["files"] = {"H300": [H300.name, sha(H300)], "NS": [NS.name, sha(NS)], "EU": [EU.name, sha(EU)],
                    "V5D_M (作り直し)": ["problem_d155_euler_v5d_M.yaml", V5D_M_SHA], "CSV": [str(CSV.relative_to(HERE)), sha(CSV)]}
    if NS.read_text() != ns_text():
        fails.append(f"{NS.name} が make の生成と違う (手で編集された?)")
    if EU.read_text() != eu_text(D["NS"]):
        fails.append(f"{EU.name} が make の生成と違う (手で編集された?)")
    if sha(CSV) != CSV_SHA:
        fails.append(f"δ_r の CSV の sha256 {sha(CSV)[:16]}… が登録 {CSV_SHA[:16]}… と違う")

    def pair(x, y, allowed, want):
        d = {p: v for p, v in diff_paths(D[x], D[y]).items() if p != ("name",)}
        keys = {".".join(map(str, p)) for p in d}
        rec[f"{x} vs {y}"] = {".".join(map(str, p)): [v[0], v[1]] for p, v in d.items()}
        extra, lack = sorted(keys - set(allowed)), sorted(set(allowed) - keys)
        if extra or lack:
            fails.append(f"{x} と {y} の差が登録と違う: 余分 {extra}、足りない {lack}")
        for p, v in want.items():
            got = d.get(tuple(p.split(".")))
            if got != v:
                fails.append(f"{x} と {y} の {p} = {got} ({v} であること)")

    ns_ev = D["NS"]["evaluate"]
    pair("H300", "NS", NS_ALLOWED, {"spec.wall_thermal.Tw": (300.0, TW_NEW), "evaluate.cfl_main": (1.0, MISSING),
                                    "evaluate.implicit_relax": (MISSING, 0.7), "evaluate.nStepOuter": (12000, 200000),
                                    "evaluate.outStepInterval": (1000, 2500), "evaluate.main_line": (MISSING, ns_ev.get("main_line")),
                                    "evaluate.stop_rule": (MISSING, True), "evaluate.stop_rule_cap": (MISSING, 200000),
                                    "evaluate.stop_rule_hold_steps": (MISSING, 40000),
                                    "evaluate.stop_rule_thresholds": (MISSING, STOP_RULE_THRESHOLDS)})
    ns_ev_nosr = {k: v for k, v in ns_ev.items() if k not in STOP_RULE_KEYS}
    for ref in ("G0C_CHAIN", "C50"):
        if ns_ev_nosr != D[ref]["evaluate"]:
            fails.append(f"NS の evaluate (停止規則のキーを除く) が {ref} の evaluate と違う: {diff_paths(D[ref]['evaluate'], ns_ev_nosr)}")
    pair("V5D_M", "EU", EU_ALLOWED, {"spec.r_throat": (D["V5D_M"]["spec"]["r_throat"], D["NS"]["spec"]["r_throat"]),
                                     "geometry.physical_upstream": (MISSING, D["NS"]["geometry"]["physical_upstream"]),
                                     "mesh_euler.throat_refine": (4.0, EU_THROAT_REFINE), "evaluate.nStepOuter": (12000, EU_NSTEP)})
    eu_ev = D["EU"]["evaluate"]
    if (eu_ev.get("outStepInterval"), eu_ev.get("cfl_main")) != (1000, 2.0):
        fails.append(f"Euler の evaluate の出力の間隔・cfl_main が登録 (1000・2.0) と違う: {(eu_ev.get('outStepInterval'), eu_ev.get('cfl_main'))}")
    # Euler と NS の設計の入力 (gas・spec・dv・geometry) が同じ (NS だけのキーを除く)
    cross = {}
    for blk in ("gas", "spec", "dv", "geometry"):
        a, b = dict(D["EU"].get(blk) or {}), dict(D["NS"].get(blk) or {})
        if blk == "spec":
            b.pop("wall_thermal", None)
        if blk == "geometry":
            for k in NS_ONLY_GEOM:
                b.pop(k, None)
                if k in a:
                    fails.append(f"Euler の geometry に NS だけのキー {k} がある")
        dd = diff_paths(a, b)
        cross[blk] = {".".join(map(str, p)): [v[0], v[1]] for p, v in dd.items()}
        if dd:
            fails.append(f"Euler と NS の {blk} が違う: {cross[blk]}")
    rec["EU vs NS (設計の入力)"] = cross
    # forge_design で読む: main_line の解決・壁温・physical_upstream
    from forge_design.evaluate import runner_axismach as RA
    pn, pe = RA.load_problem(NS), RA.load_problem(EU)
    ml = RA._main_line_spec(pn)
    rec["main_line_resolved"] = ml
    want_ml = {"line": "dir", "cap": 50.0, "itj": 5, "cfl": 4.0, "relax": 0.7, "prestage_relax": 0.7,
               "limiter_ref": {"limiterRoRef": 0.8739869154, "limiterPRef": 387491.0573, "limiterARef": 359.7712185},
               "extra_fields": ["res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega", "omg_prod", "omg_dest", "omg_cross",
                                "omg_trans", "omg_axisym", "dt_local"]}
    if ml != want_ml:
        fails.append(f"main_line の解決が登録と違う: {ml}")
    if pn.wall_thermal != {"mode": "isothermal", "Tw": TW_NEW}:
        fails.append(f"NS の wall_thermal が 370 K の等温でない: {pn.wall_thermal}")
    rec["wall_thermal"] = pn.wall_thermal
    rec["physical_upstream"] = {"ns": RA._physical_upstream(pn.geometry), "euler": RA._physical_upstream(pe.geometry)}
    if rec["physical_upstream"]["ns"] != rec["physical_upstream"]["euler"]:
        fails.append("NS と Euler の physical_upstream の解決が違う")
    rec["mesh_euler"] = RA.mesh_euler_block(pe)
    if (rec["mesh_euler"].get("ni"), rec["mesh_euler"].get("nj"), rec["mesh_euler"].get("throat_refine")) != (2000, 97, EU_THROAT_REFINE):
        fails.append(f"Euler の mesh_euler が 2000 × 97・throat_refine {EU_THROAT_REFINE} でない: {rec['mesh_euler']}")
    rec["euler_evaluate"] = {k: eu_ev.get(k) for k in ("nStepOuter", "outStepInterval", "cfl_main", "mid_stage")}
    # 停止規則のキーの解決 (stop_rule.resolve_settings; euler_ref は投入の時に固定するので、ここでは形だけの仮の値)
    from forge_design.evaluate import stop_rule as SR
    n, out_int, n_src = RA._ns_steps(pn, None)
    sr = SR.resolve_settings(pn.evaluate, nsteps_arg=None, n_main=n, n_main_source=n_src, out_int=out_int,
                             euler_ref="<投入の時に固定する run_0950 の res と sha256>", requested=None)
    rec["stop_rule_resolved"] = sr
    want_sr = {"use": True, "requested": True, "cap": 200000, "hold_steps": 40000, "out_step_interval": 2500, "fallback_steps": 200000,
               **{k: (list(v) if isinstance(v, list) else v) for k, v in STOP_RULE_THRESHOLDS.items()}}
    bad_sr = {k: sr.get(k) for k, v in want_sr.items() if sr.get(k) != v}
    if bad_sr:
        fails.append(f"停止規則の解決が登録と違う: {bad_sr}")
    rec["fails"] = fails
    rec["VERDICT"] = "PASS" if not fails else "FAIL"
    OUT.mkdir(parents=True, exist_ok=True)
    (Path(out_json) if out_json else OUT / "h370_problems.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    print("\n".join(fails) if fails else "差は登録のとおり")
    print(f"VERDICT: {rec['VERDICT']}")
    return 0 if not fails else 1


# --- prep-check (CFD 0 step) ------------------------------------------------------------------------------------------------
def load_pinned(path: Path, pin: Path, td: Path):
    from forge_design.evaluate import runner_axismach as RA
    raw = yaml.safe_load(path.read_text())
    raw["geometry"]["initial_line_run"] = str(pin)
    q = td / f"_pin_{path.name}"
    q.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    return RA.load_problem(q)


def generator_sha256(coords, quads, bedges) -> str:
    """runner_axismach._mesh_hashes の generator_sha256 と同じ式 (変換器に依らない生成器の出力の sha256)。"""
    h = hashlib.sha256()

    def _upd(name, a):
        a = np.ascontiguousarray(a)
        h.update(name.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    _upd("coords", np.asarray(coords, dtype=np.float64))
    _upd("quads", np.asarray(quads, dtype=np.int64))
    for g in ("inlet", "outlet", "wall", "axis"):
        _upd(g, np.asarray(bedges[g], dtype=np.int64))
    return h.hexdigest()


def quality(coords, quads, ar_max: float) -> dict:
    from check_mesh_quality import metrics_vectorized
    nq = len(quads)
    ars, sk = metrics_vectorized(coords[:, :2], quads.ravel(), np.arange(1, nq + 1) * 4, np.full(nq, 9), False)
    i = int(np.argmax(ars))
    c = coords[quads[i], :2].mean(axis=0)
    return {"cells": int(nq), "ar_max": float(ars.max()), "ar_p99": float(np.percentile(ars, 99)), "skew_max": float(sk.max()),
            "n_ar_gt_max": int((ars > ar_max).sum()), "n_skew_gt_0.9": int((sk > 0.9).sum()), "ar_limit": float(ar_max),
            "at_ar_max_xy_m": [float(c[0]), float(c[1])], "pass": bool(ars.max() <= ar_max and sk.max() <= 0.9),
            "note": "生成の座標に当てた値 (変換後の HDF5 の check_mesh_quality の VERDICT ではない — AWS で確かめる)"}


def first_layer(coords, ni: int, scale: float) -> dict:
    P = coords[:, :2].reshape(ni, -1, 2)
    d = np.hypot(*(P[:, -1] - P[:, -2]).T)
    x = P[:, -1, 0] / scale
    out = {"min_m": float(d.min()), "x_at_min": float(x[int(np.argmin(d))]), "max_m": float(d.max()), "x_at_max": float(x[int(np.argmax(d))])}
    for lab, (lo, hi) in {"x_lt_0": (-np.inf, 0.0), "throat_0_1": (0.0, 1.0), "test_40_94": (40.0, 94.0)}.items():
        m = (x >= lo) & (x < hi)
        out[lab] = [float(d[m].min()), float(d[m].max())]
    out["note"] = "壁の節点と第一の内部の節点の距離 (壁法線の近壁層では法線方向)。変換器の wall_dist (双対の重心の距離) は AWS で"
    return out


class HRefWall:
    """診断の参照の壁 (Euler の参照の壁の定義; #34 の A と同じ): 直管 r_inlet ([x_in, −L_U))・(−L_U, r_inlet, 0, 0) → x = 0 の S の
    (r, r′, r″) の 5 次 Hermite ([−L_U, 0))・S (x ≥ 0)。#36e からは正式の経路 (euler_reference_wall の案 A) と一致することの確認に使う。"""

    def __init__(self, design, r_U, L_U):
        from forge_design.geometry.wall_walldriven import _hermite_quintic
        self.d, self.r_U, self.L_U = design, float(r_U), float(L_U)
        self.x_in, self.x_e = float(design.x_in), float(design.x_e)
        self.c = _hermite_quintic(-self.L_U, 0.0, self.r_U, 0.0, 0.0, *[float(design.r(np.r_[0.0], n)[0]) for n in range(3)])

    def r(self, x, deriv: int = 0):
        from forge_design.geometry.wall_walldriven import _poly_eval
        x = np.asarray(x, dtype=float)
        out = np.array(self.d.r(x, deriv), dtype=float, copy=True)
        mp, mq = x < -self.L_U, (x >= -self.L_U) & (x < 0.0)
        out[mp] = self.r_U if deriv == 0 else 0.0
        if mq.any():
            out[mq] = _poly_eval(self.c, (x[mq] + self.L_U) / self.L_U, deriv) / self.L_U ** deriv
        return out


def cmd_prep_check(pin: Path) -> int:
    from forge_design.evaluate import runner_axismach as RA
    from forge_design.geometry.wall_axismach import check_required_attrs, monotone_gate
    from forge_design.meshing.mesh2d import generate_axisym_mesh
    from forge_design.metrics.delta_in_check import compare_delta_in
    fails = []
    rec = {"plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36b・#36c・#36e", "tool": "h370_problems.py prep-check",
           "environment": RA.environment_record(), "cfd_pin_dir": str(pin), "files": {NS.name: sha(NS), EU.name: sha(EU), "csv": sha(CSV)},
           "scope": "CFD 0 step: prepare_ns の変換の前までを同じ関数で再現 (変換器・forge は使わない)"}
    if sha(CSV) != CSV_SHA:
        raise SystemExit("δ_r の CSV の sha256 が登録と違う — 止める")
    td = Path(tempfile.mkdtemp(prefix="h370_"))
    # --- NS (prepare_ns と同じ順) ---
    pn = load_pinned(NS, pin, td)
    main_line = RA._main_line_spec(pn, None, None)
    pu = RA._physical_upstream(pn.geometry)
    pwu = RA._pw_upstream(pn.geometry)
    d = RA.design_chain(pn)
    RA.require_moc_gate(pn, d)
    scale = float(pn.spec["r_throat"])
    tbl = np.loadtxt(CSV, delimiter=",", skiprows=1)
    drx = RA.delta_r_from_table(tbl[:, 0], tbl[:, 1])
    wall = RA.build_physical_wall(pn, d, scale, delta_r_x=drx, offset="radial", pwu=pwu)
    check_required_attrs(wall)
    msgs = wall.validate()
    if msgs:
        fails.append(f"NS の物理壁の validate: {msgs}")
    rec["design"] = {"cfd_pin_match": (d.get("initial_line") or {}).get("match"), "x_in": float(d["wall"].x_in), "x_e": float(d["wall"].x_e),
                     "n_wall_inv": int(len(d["wall_inv"])), "L_c": float(d["L_c"]), "moc_gate": (d.get("moc") or {}).get("gate", {}).get("pass")}
    rec["ns_wall"] = {"type": type(wall).__name__, "physical_upstream": wall.physical_upstream,
                      "throat": {"x": float(wall.x_throat), "r": float(wall.r_throat), "kappa": float(wall.kappa_throat)},
                      "upstream_gate_pass": wall.upstream_gate.get("pass")}
    mp = RA.mesh_params(pn, scale, ni=561, nj=97, wall_first_frac=4.5e-5)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    ni = mp.ni
    nj = coords.shape[0] // ni
    rec["ns_mesh"] = {"ni": ni, "nj": nj, "axis_cap_frac": mp.axis_cap_frac, "wall_normal_layer": mp.wall_normal_layer,
                      "msh_digits": int(pn.mesh.get("msh_digits", 10)), "generator_sha256": generator_sha256(coords, quads, bedges),
                      "quality": quality(coords, quads, float(pn.mesh.get("ar_max", 1000))), "first_layer": first_layer(coords, ni, scale)}
    if not rec["ns_mesh"]["quality"]["pass"]:
        fails.append(f"NS の格子の品質 (生成の座標) が不合格: {rec['ns_mesh']['quality']}")
    if (ni, nj, mp.axis_cap_frac) != (4719, 160, 0.015):
        fails.append(f"NS の格子が Gref (4719 × 160・axis_cap_frac 0.015) でない: {(ni, nj, mp.axis_cap_frac)}")
    # 一定 x の断面への補間の端の補間範囲の不足 (受け渡しの抽出の remap_edge="missing" で欠測になる列; 0 であること。#36e (4))
    Xg = coords[:, 0].reshape(ni, nj) / scale
    xt_ = Xg[:, 0]
    miss_cols = np.zeros(ni, dtype=bool)
    for j in range(nj):
        miss_cols |= (xt_ < Xg[0, j]) | (xt_ > Xg[-1, j])
    rec["ns_mesh"]["remap_constant_x"] = {"n_missing_cols": int(miss_cols.sum()), "x_missing": [float(v) for v in xt_[miss_cols][:20]],
                                          "max_column_dx_rt": float(np.max(np.ptp(Xg, axis=1))),
                                          "note": "生成の座標で、行 j の x の範囲が列の目標 x[:, 0] を覆わない列 (抽出で欠測になる)"}
    if miss_cols.any():
        fails.append(f"NS の格子で一定 x の断面への補間の範囲が足りない列 {int(miss_cols.sum())} (受け渡しの抽出で欠測になる)")
    # 停止規則の要求の断面の位置 (x_fracs × x_w、x_w = 物理壁の出口) と最も近い列 (§6 V-h370)
    x_w = float(wall.x_e)
    rec["stop_rule_sections"] = [{"x_frac": f, "x_requested": f * x_w, "x_column": float(xt_[int(np.argmin(np.abs(xt_ - f * x_w)))])}
                                 for f in STOP_RULE_X_FRACS]
    # 300 K の H (#34) と格子・壁がビット同一 (370 K の問題は壁温と evaluate だけが違う)
    p3 = load_pinned(H300, pin, td)
    d3 = RA.design_chain(p3)
    w3 = RA.build_physical_wall(p3, d3, scale, delta_r_x=drx, offset="radial")
    c3, q3, b3 = generate_axisym_mesh(w3, RA.mesh_params(p3, scale, ni=561, nj=97, wall_first_frac=4.5e-5))
    rec["ns_mesh"]["same_as_tw300_h"] = bool(np.array_equal(c3, coords) and np.array_equal(q3, quads))
    if not rec["ns_mesh"]["same_as_tw300_h"]:
        fails.append("NS の格子が 300 K の H (#34) とビット同一でない")
    # 設定 (prepare_ns と同じ引数の解決)
    n = int(pn.evaluate.get("nStepOuter", 48000))
    out_int = int(pn.evaluate.get("outStepInterval", max(n // 6, 1)))
    cfl_main = float(pn.evaluate.get("cfl_main", 1.0))
    implicit_relax = pn.evaluate.get("implicit_relax")
    bc = RA._bcond_with_species(RA._bcond(pn, euler=False, profile_ok=True), RA._tp_species_Y(pn))
    wall_lines = [ln for ln in bc.splitlines() if "wall" in ln]
    rec["bcond_wall_lines"] = wall_lines
    if not any("wall_isothermal" in ln and "370.0" in ln for ln in wall_lines):
        fails.append(f"bcond の壁の行に wall_isothermal・370.0 が無い: {wall_lines}")
    rd = td / "cfg"
    rd.mkdir()
    cfg_ns, cfg_main = RA._ns_configs(pn, rd, n, out_int, cfl_main, implicit_relax, main_line)
    pre, cfg_main_eff, main_want = RA.ns_stage_configs(cfg_main, "full", cfg_pre=cfg_ns)
    keys = {"cfl": ("time", "deltaT", "cfl"), "cfl_pseudo": ("time", "deltaT", "cfl_pseudo"), "convMethod": ("space", "convMethod"),
            "nStepOuter": ("time", "last", "nStepOuter"), "outStepInterval": ("time", "outStepInterval"),
            "implicitRelax": ("time", "deltaT", "implicitRelax"), "lineImplicit": ("time", "deltaT", "lineImplicit"),
            "lineDtDirectional": ("time", "deltaT", "lineDtDirectional"), "lineDtDirectionalCap": ("time", "deltaT", "lineDtDirectionalCap"),
            "implicitThermalJacobian": ("time", "deltaT", "implicitThermalJacobian"), "limiterRoRef": ("space", "limiterRoRef"),
            "limiterPRef": ("space", "limiterPRef"), "limiterARef": ("space", "limiterARef"), "output": ("output",),
            "blockDPLUR": ("time", "deltaT", "blockDPLUR"), "timeIntegration": ("time", "timeIntegration"),
            "nStepInner": ("time", "nStepInner"), "limiter": ("space", "limiter")}

    def eff(cfg):
        doc = RA._cfg_load(cfg)
        return {k: (lambda v: None if v is RA._CFG_MISSING else v)(RA._cfg_path_get(doc, pth)) for k, pth in keys.items()}
    stages = {tag: eff(cfg) for tag, cfg, _, _ in pre}
    stages["main"] = eff(cfg_main_eff)
    rec["stages"] = stages
    rec["stage_sha256"] = RA.ns_stage_sha256(cfg_main, cfg_ns)
    want = {"S1_soft": {"cfl": 0.5, "convMethod": 0, "nStepOuter": 3000}, "S2_mid": {"cfl": 1.0, "convMethod": 0, "nStepOuter": 3000},
            "main": {"cfl": 4.0, "nStepOuter": 200000, "outStepInterval": 2500, "lineImplicit": 1, "lineDtDirectional": 1,
                     "lineDtDirectionalCap": 50.0, "implicitThermalJacobian": 5, "limiterRoRef": 0.8739869154, "limiterPRef": 387491.0573,
                     "limiterARef": 359.7712185}}
    for tag, w in want.items():
        for k, v in w.items():
            if stages.get(tag, {}).get(k) != v:
                fails.append(f"{tag} の {k} = {stages.get(tag, {}).get(k)!r} ({v!r} であること)")
    for tag in ("S1_soft", "S2_mid"):
        left = [k for k in ("lineImplicit", "lineDtDirectional", "lineDtDirectionalCap", "implicitThermalJacobian", "limiterRoRef",
                            "limiterPRef", "limiterARef") if stages[tag][k] is not None]
        if left or stages[tag]["output"] is not None:
            fails.append(f"{tag} に B0 のキー {left} (output {stages[tag]['output']}) が残っている")
    if not all(stages[t]["implicitRelax"] == 0.7 for t in stages):
        fails.append(f"implicitRelax が 3 段で 0.7 でない: {[stages[t]['implicitRelax'] for t in stages]}")
    # --- Euler の参照 ---
    pe = load_pinned(EU, pin, td)
    de = RA.design_chain(pe)
    RA.require_moc_gate(pe, de)
    S, Sn = de["wall"], d["wall"]
    scale_e = float(pe.spec["r_throat"])
    xs = np.r_[0.0, np.linspace(0.0, float(S.x_e), 400001)[1:]]
    same_S = (float(S.x_e) == float(Sn.x_e) and np.array_equal(de["wall_inv"], d["wall_inv"])
              and all(np.array_equal(S.r(xs, k), Sn.r(xs, k)) for k in range(4)))
    rec["euler_design"] = {"same_S_as_ns": bool(same_S), "scale_same": bool(scale_e == scale), "cfd_pin_match": (de.get("initial_line") or {}).get("match"),
                           "monotone_gate_S": {k: monotone_gate(S._spl, S.x0, S.x_e)[k] for k in ("pass", "min_r1", "x_min_theta", "D", "x_D_peak", "x_D")}}
    if not same_S:
        fails.append("Euler と NS の設計壁 S がビット同一でない")
    # 正式の Euler の参照の壁 (#36e の案 A: H の上流を付けた設計壁)。ゲートの不合格は例外
    from forge_design.geometry.wall_axismach import DESIGN_UPSTREAM_KIND, load_wall_file, save_wall_file
    wE, er = RA.euler_reference_wall(pe, de, scale_e)
    rec["euler_reference_wall"] = {"status": "built", "type": type(wE).__name__, "repr": er["repr"], "throat": er["throat"],
                                   "gate": {k: er["gate"][k] for k in ("pass", "reasons", "downstream_monotone_gate", "upstream", "max_seam_jump",
                                                                       "max_abs_d2_change_vs_H")},
                                   "global_min": {k: er["gate"]["global_min"][k] for k in ("pass", "x_t", "r_t", "unique", "monotone_after")}}
    if not er["gate"]["pass"]:
        fails.append(f"Euler の参照の壁のゲートが不合格: {er['gate']['reasons']}")
    Wd = HRefWall(S, pu["r_inlet"], pu["L_U"])           # 診断の参照の壁 (#34 の A の定義) との一致
    xa = np.linspace(float(S.x_in), -1e-13, 250001)
    rec["euler_reference_wall"]["vs_HRefWall"] = {"x_ge_0_bit_identical": bool(all(np.array_equal(wE.r(xs, k), Wd.r(xs, k)) for k in range(3))),
                                                  "x_lt_0_max_abs": float(max(np.max(np.abs(wE.r(xa, k) - Wd.r(xa, k))) for k in range(3)))}
    mpE = RA.mesh_params_euler(pe, scale_e)
    cE, qE, bE = generate_axisym_mesh(wE, mpE)
    niE = mpE.ni
    mE = RA.mesh_euler_block(pe)
    rec["euler_mesh"] = {"wall": "正式の Euler の参照の壁 (euler_reference_wall、案 A)", "ni": niE, "nj": cE.shape[0] // niE, "mesh_euler": mE,
                         "throat_refine": mpE.throat_refine, "generator_sha256": generator_sha256(cE, qE, bE),
                         "quality": quality(cE, qE, float(mE.get("ar_max", 1000)))}
    if not rec["euler_mesh"]["quality"]["pass"]:
        fails.append(f"Euler の格子の品質 (生成の座標) が不合格: {rec['euler_mesh']['quality']}")
    if (niE, rec["euler_mesh"]["nj"], mpE.throat_refine) != (2000, 97, EU_THROAT_REFINE):
        fails.append(f"Euler の格子が 2000 × 97・throat_refine {EU_THROAT_REFINE} でない: {(niE, rec['euler_mesh']['nj'], mpE.throat_refine)}")
    # 壁ファイルの保存 → 復元 → 同じ格子の生成: 座標・接続がビット一致 (#36e (1))
    wdir = td / "wallfile"
    wdir.mkdir()
    wpath, _, wsha = save_wall_file(wdir, wE, scale_e, DESIGN_UPSTREAM_KIND)
    back = load_wall_file(wpath)
    cB, qB, bB = generate_axisym_mesh(back["physical"], mpE)
    rec["euler_wall_file"] = {"kind": back["kind"], "sha256": wsha, "coords_bit_identical": bool(np.array_equal(cB, cE)),
                              "quads_identical": bool(np.array_equal(qB, qE)),
                              "bedges_identical": bool(all(np.array_equal(np.asarray(bB[g]), np.asarray(bE[g])) for g in bE))}
    if not (rec["euler_wall_file"]["coords_bit_identical"] and rec["euler_wall_file"]["quads_identical"]
            and rec["euler_wall_file"]["bedges_identical"]):
        fails.append(f"壁ファイルから復元した壁の格子が元とビット一致しない: {rec['euler_wall_file']}")
    XN = coords[:, 0].reshape(ni, nj) / scale
    RN = coords[:, 1].reshape(ni, nj) / scale
    XE = cE[:, 0].reshape(niE, -1) / scale_e
    RE = cE[:, 1].reshape(niE, -1) / scale_e
    res = compare_delta_in(XN[:, 0], RN[:, -1], XE[:, 0], RE[:, -1], drx, scale, r_ns_exact=wall.r, r_eu_exact=wE.r, return_arrays=True)
    arr = res.pop("arrays")
    rec["delta_in_check"] = res
    # 列の配置と補間の誤差の位置: x ≥ 0 で tol を超える列の範囲・数
    err = arr["delta_in"] - arr["delta_r"]
    over = (arr["x"] >= 0.0) & (np.abs(err) * scale * 1e3 > res["tol_mm"])
    rec["delta_in_check"]["over_tol"] = {"n": int(over.sum()), "x_range": ([float(arr["x"][over].min()), float(arr["x"][over].max())]
                                                                         if over.any() else None)}
    np.savetxt(OUT / "h370_delta_in_tr5.csv", np.c_[arr["x"], arr["delta_in"], arr["delta_r"], err * scale * 1e3], delimiter=",",
               header="x_rt,delta_in_rt,delta_r_rt,err_mm", comments="")
    # 旧い参照 (run_0174、出口の終点の修正の前) の終端 95.2275815792169 の外の NS の列 (A/B で欠ける区間)
    x_e_old = 95.2275815792169
    rec["old_reference_uncovered_ns_cols"] = {"x_e_old": x_e_old, "n": int(np.sum(XN[:, 0] > x_e_old + 1e-9)),
                                              "x_range": [float(XN[XN[:, 0] > x_e_old + 1e-9, 0].min()), float(XN[-1, 0])]}
    if not res["pass"]:
        fails.append(f"δ_in の照合が不合格 (x ≥ 0 の max {res['judged'].get('max_abs_mm')} mm > {res['tol_mm']} mm、被覆 {res['covered']})")
    rec["fails"] = fails
    rec["VERDICT"] = "PASS" if not fails else "FAIL"
    rec["aws_todo"] = ["prepare_ns (変換器 2 回) の後の MESH_QUALITY.txt (check_mesh_quality --ar-max 5000) の VERDICT と高 AR のセルの位置",
                       "変換後の nozzle.h5 の /VALUE/wall_dist (双対の重心の距離) と第一層の局所の法線の距離・y₁⁺ の見込み",
                       "prepare_info.json の mesh.hashes.generator_sha256 がここの値と一致 (同じ環境なら; 違えば numpy・scipy の版の差を記録)",
                       "solverConfig_prestage.yaml・solverConfig.yaml・ns_stage_configs の sha256 とここの stage_sha256 の一致",
                       "Euler の参照: prepare (変換器・品質・IC) の MESH_QUALITY.txt の VERDICT と、check_runs (delta_in_check の CLI) で"
                       "変換後の座標の δ_in の照合 (≤ 0.001 mm)",
                       "Euler の参照の資格を満たした時点の res と sha256 を固定して NS に渡す (deltastar_loop --euler-ref-res・--euler-ref-sha256)"]
    rec["superseded_record"] = {"file": "h370_prep_check.json・h370_delta_in.csv",
                                "note": "throat_refine 4・診断の参照の壁 (HRefWall) での照合の旧い FAIL (max 1.1367e-3 mm)。消さずに残す (§5.1 #36e)"}
    (OUT / "h370_prep_check_tr5.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    print(json.dumps({k: rec[k] for k in ("design", "ns_wall", "euler_design", "euler_reference_wall", "euler_wall_file",
                                          "old_reference_uncovered_ns_cols", "stop_rule_sections")},
                     ensure_ascii=False, indent=1, default=str)[:6000])
    print(json.dumps({"ns_mesh": {k: rec["ns_mesh"][k] for k in ("ni", "nj", "quality", "first_layer", "same_as_tw300_h", "generator_sha256",
                                                                "remap_constant_x")},
                      "euler_mesh_quality": rec["euler_mesh"]["quality"], "stages": rec["stages"],
                      "delta_in": {k: rec["delta_in_check"][k] for k in ("pass", "covered", "judged", "by_region", "decomposition", "over_tol")}},
                     ensure_ascii=False, indent=1, default=str))
    print("\n".join(fails) if fails else "")
    print(f"VERDICT: {rec['VERDICT']}")
    return 0 if not fails else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("make", "check", "prep-check"))
    ap.add_argument("--pin", default=None)
    a = ap.parse_args()
    if a.cmd == "make":
        return cmd_make()
    if a.cmd == "check":
        return cmd_check()
    pin = Path(a.pin).resolve() if a.pin else PIN_DEFAULT
    if not pin.is_dir():
        raise SystemExit(f"CFD ピンの凍結元 {pin} が無い — --pin で渡す")
    return cmd_prep_check(pin)


if __name__ == "__main__":
    sys.exit(main())
