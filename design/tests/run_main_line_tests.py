#!/usr/bin/env python3
"""`evaluate.main_line` (本段だけ B0) と NS の段の設定の導出の試験 (plan tooling-nozzle-isothermal-wall-chain §4.8-2・§5.1 #30)。

(a) case/45 の tw300 の YAML に main_line (cap 50、リミッタの基準値は run_0752・run_0717 と同じ、extra_fields は run_0752 と同じ 12 量) を
    書いた写しで作る本段の solverConfig が、run_0752_t2_d1/solverConfig.yaml (AWS の実物、上限だけ 1.0。下の REF_0752) と
    nStepOuter・outStepInterval・lineDtDirectionalCap・mesh.axisSegmentRWeight (相手は既定と同じ 1 を明示) 以外でキーごとに一致する。
    その 4 つを揃えると文字列でも一致する (書き方が cold_cfl.py prep と同じ)。
(b) S1_soft・S2_mid の設定は B0 を当てる前の設定から作られ、ライン陰解法・上限・ITJ・リミッタの基準値・extraFields が無く、
    cfl・次数・nStepInner は今と同じ (main_line の無い元の YAML の S1・S2 と文字列で一致)。
(c) main_line の無い問題は前の段の元の設定を書かず、段の設定の導出は従来と同じ (固定の入力に対する文字列の期待値)。
    変更前のコード (HEAD 09d1d39f) との文字列一致は git worktree で別に確かめた (§5.1 #30 の行)。
(d) run_staged_ns を forge の代わりに差し替えて回し、S1・S2・本段が順に書かれ、通知 MAIN_STAGE_STARTED が前の段の出力の処理の後・
    本段の起動の前に出て、本段の設定・初期の場・前の段の設定の sha256 が実際のファイルと一致する。
(e) 不正な main_line・cfl_main との併記・main_line.relax と違う implicit_relax・ramp との組み合わせ・前の段の元の設定と本段の不整合は ValueError。
design_chain・メッシュ・変換器・forge は使わない (設定を作る関数と段の導出を関数の単位で試す)。
"""
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0
SKIPPED = []
TW300 = ROOT / "case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300.yaml"
TMP = Path(tempfile.mkdtemp(prefix="main_line_"))
ML = {"line": "dir", "cap": 50, "itj": 5, "cfl": 4, "relax": 0.7,
      "limiter_ref": {"ro": 0.8739869154, "p": 387491.0573, "a": 359.7712185},
      "extra_fields": ["res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega",
                       "omg_prod", "omg_dest", "omg_cross", "omg_trans", "omg_axisym", "dt_local"]}
# AWS の case/45.isobutane_m6_d155/run_0752_t2_d1/solverConfig.yaml の写し (sha256 4e4fa4331feae32f…; cold_cfl.py prep
# --cfl 4 --limiter-ref-from … --line dir --itj 5 --cap 1 --extra … + mesh.axisSegmentRWeight 1 の 1 行)
REF_0752 = '''mesh: {axisSegmentRWeight: 1, discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 2, viscMethod: 2,
           visc: 1.8e-5, thermCond: 0.0257, prandtlLam: 0.72, cp: 1360.0, gamma: 1.27354,
           species: [{name: "MIXDRY", lump: {"CO2": 0.12550903682310086, "O2": 0.022645262404230745, "N2": 0.8518457007726683}, basis: mole}, "H2O"],
           transport: {"CO2": "cea", "O2": "cea", "N2": "cea", "H2O": "custom:h2o_iapws_cea_v1"}, thermoHrefTemp: 298.15}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 170000}
  deltaT: {lineDtDirectionalCap: 1.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1, control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0,
           dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 2500
  timeIntegration: 11
  nStepInner: 5
space: {limiterRoRef: 0.8739869154, limiterPRef: 387491.0573, limiterARef: 359.7712185, convMethod: 1, limiter: 2}
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0, turbulentPrandtl: 0.9}
initial: "uniform_p101325_u10"
output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}
'''
B0_ONLY = [("time", "deltaT", k) for k in ("lineImplicit", "lineDtDirectional", "lineDtDirectionalCap", "implicitThermalJacobian")] \
    + [("space", k) for k in ("limiterRoRef", "limiterPRef", "limiterARef")] + [("output",)]


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


def rejects(name, fn, must=""):
    try:
        fn()
    except ValueError as e:
        check(f"{name} を拒否", must in str(e), str(e)[:110])
        return
    check(f"{name} を拒否 (通ってしまった)", False)


def problem(tag, edit):
    """tw300 の YAML の写しに edit(raw) を当てて書き、load_problem で読む。"""
    raw = yaml.safe_load(TW300.read_text())
    edit(raw)
    path = TMP / f"{tag}.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    return load_problem(path)


def with_ml(ml=None, drop_cfl=True, **ev):
    def edit(raw):
        if drop_cfl:
            raw["evaluate"].pop("cfl_main", None)
        if ml is not None:
            raw["evaluate"]["main_line"] = ml
        raw["evaluate"].update(ev)
    return edit


def configs(p, tag, cfl_main=None, implicit_relax=None, ml=None):
    rd = TMP / f"rd_{tag}"
    rd.mkdir()
    n = int(p.evaluate.get("nStepOuter", 48000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    return RA._ns_configs(p, rd, n, out_int, cfl_main, implicit_relax, ml)


def get(text, path):
    return RA._cfg_path_get(RA._cfg_load(text), path)


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


if not TW300.is_file():
    SKIPPED.append(f"(a)〜(e) ({TW300.relative_to(ROOT)} が無い)")
else:
    # --- (a) 本段の設定 = run_0752 の実物 (許した 4 項目以外) --------------------------------------------------------
    p_ml = problem("ml", with_ml(ML))
    ml = RA._main_line_spec(p_ml)
    check("(a) main_line を解決 (cap 50.0・itj 5・cfl 4.0・relax 0.7・基準値 3 つ・出力 12 量)",
          ml["cap"] == 50.0 and ml["itj"] == 5 and ml["cfl"] == 4.0 and ml["relax"] == 0.7
          and ml["limiter_ref"] == {"limiterRoRef": 0.8739869154, "limiterPRef": 387491.0573, "limiterARef": 359.7712185}
          and len(ml["extra_fields"]) == 12, ml)
    base, main = configs(p_ml, "ml", ml=ml)
    ref = REF_0752
    diff = RA._cfg_diff_paths(RA._cfg_load(ref), RA._cfg_load(main))
    allowed = {("time", "last", "nStepOuter"), ("time", "outStepInterval"), ("time", "deltaT", "lineDtDirectionalCap"),
               ("mesh", "axisSegmentRWeight")}
    check("(a) run_0752 との差は nStepOuter・outStepInterval・lineDtDirectionalCap・mesh.axisSegmentRWeight だけ",
          set(diff) <= allowed, {".".join(k): v for k, v in diff.items() if k not in allowed})
    check("(a) その 4 つの値 (相手 170000・2500・1.0・1 / こちら 12000・1000・50.0・無し = 既定 1)",
          diff.get(("time", "last", "nStepOuter")) == (170000, 12000) and diff.get(("time", "outStepInterval")) == (2500, 1000)
          and diff.get(("time", "deltaT", "lineDtDirectionalCap")) == (1.0, 50.0)
          and diff.get(("mesh", "axisSegmentRWeight"), (None, None))[0] == 1
          and diff.get(("mesh", "axisSegmentRWeight"), (None, None))[1] is RA._CFG_MISSING, diff)
    norm = (main.replace("nStepOuter: 12000", "nStepOuter: 170000").replace("outStepInterval: 1000", "outStepInterval: 2500")
            .replace("lineDtDirectionalCap: 50.0", "lineDtDirectionalCap: 1.0").replace("mesh: {", "mesh: {axisSegmentRWeight: 1, ", 1))
    check("(a) 4 つを揃えると run_0752 と文字列でも一致 (キーの順・値の書き方が cold_cfl.py prep と同じ)", norm == ref)
    # 前の段の元の設定 = main_line の無い元の YAML の本段の設定 (cfl_main 1.0 は既定と同じなので外しても同じ)
    p_orig = load_problem(TW300)
    base0, main0 = configs(p_orig, "orig")
    check("(a) 前の段の元の設定 = main_line の無い元の YAML の設定 (文字列一致)", base == base0 == main0)
    check("(a) 本段と前の段の元の設定の差は B0 のパスだけ",
          set(RA._cfg_diff_paths(RA._cfg_load(base), RA._cfg_load(main))) <= RA._MAIN_LINE_PATHS)

    # --- (b) S1・S2 は B0 を当てる前の設定から ----------------------------------------------------------------------
    pre, main_f, main_want = RA.ns_stage_configs(main, "full", cfg_pre=base)
    pre0, main0_f, _ = RA.ns_stage_configs(main0, "full")
    check("(b) 段は S1_soft・S2_mid・本段", [t for t, *_ in pre] == ["S1_soft", "S2_mid"] and main_f == main)
    for (tag, cfg, nsteps, want), (_, cfg0, _, _) in zip(pre, pre0):
        left = [".".join(pth) for pth in B0_ONLY if get(cfg, pth) is not RA._CFG_MISSING]
        check(f"(b) {tag}: ライン陰解法・上限・ITJ・リミッタの基準値・extraFields が無い", not left, left)
        check(f"(b) {tag}: main_line の無い元の YAML の {tag} と文字列で一致 (cfl・次数・nStepInner・implicitRelax 無し)",
              cfg == cfg0)
        check(f"(b) {tag}: cfl {want['cfl']}・convMethod 0・nStepInner 10・3000 step",
              get(cfg, ("time", "deltaT", "cfl")) == want["cfl"] == get(cfg, ("time", "deltaT", "cfl_pseudo"))
              and get(cfg, ("space", "convMethod")) == 0 and get(cfg, ("time", "nStepInner")) == 10
              and get(cfg, ("time", "last", "nStepOuter")) == 3000 == nsteps
              and get(cfg, ("time", "deltaT", "implicitRelax")) is RA._CFG_MISSING)
    pre_bad, *_ = RA.ns_stage_configs(main, "full")
    check("(b) 対照: 前の段の元の設定を渡さないと S1 に B0 が残る (codex の指摘の再現)",
          get(pre_bad[0][1], ("time", "deltaT", "lineImplicit")) == 1
          and get(pre_bad[0][1], ("time", "deltaT", "lineDtDirectionalCap")) == 50.0)
    pre_n, main_n, _ = RA.ns_stage_configs(main, "none", cfg_pre=base)
    check("(b) stages none は本段だけ (B0 のまま)", pre_n == [] and main_n == main)

    # --- (c) main_line の無い問題: 前の段の元の設定は本段と同じ、段の導出は従来どおり --------------------------------
    check("(c) main_line の無い問題の _main_line_spec は None", RA._main_line_spec(p_orig) is None
          and RA._main_line_spec(p_orig, cfl_main=5.0, implicit_relax=0.7) is None)
    s1 = (main0.replace("nStepOuter: 12000", "nStepOuter: 3000").replace("cfl: 1.0, cfl_pseudo: 1.0", "cfl: 0.5, cfl_pseudo: 0.5")
          .replace("outStepInterval: 1000", "outStepInterval: 3000").replace("nStepInner: 5", "nStepInner: 10")
          .replace("convMethod: 1", "convMethod: 0"))
    s2 = s1.replace("cfl: 0.5, cfl_pseudo: 0.5", "cfl: 1.0, cfl_pseudo: 1.0")
    check("(c) S1・S2・本段 = 文字列の置換で作った期待値 (従来の導出)",
          [c for _, c, *_ in pre0] == [s1, s2] and main0_f == main0)
    b07, m07 = configs(p_orig, "orig_relax", cfl_main=5.0, implicit_relax=0.7)
    check("(c) cfl_main 5・implicit_relax 0.7 (cold_pair の呼び方) も従来どおり 1 つの設定",
          b07 == m07 and "cfl: 5.0, cfl_pseudo: 5.0" in m07 and "blockDPLUR: 1, implicitRelax: 0.7," in m07)
    pre_r, main_r, want_r = RA.ns_stage_configs(main0, "ramp")
    check("(c) ramp は従来どおり 2 次のまま cfl 1・2・3.5 と本段 9000 step",
          [t for t, *_ in pre_r] == ["R1_cfl1", "R2_cfl2", "R3_cfl3.5"] and want_r["nStepOuter"] == 9000
          and all(get(c, ("space", "convMethod")) == 1 for _, c, *_ in pre_r))

    # --- (d) run_staged_ns (forge を差し替え): 段の順・通知の時点と中身 -----------------------------------------------
    seen = []
    orig = (RA.run_forge, RA.stage_gate, RA._restart_same_mesh)

    def fake_forge(run_dir):
        run_dir = Path(run_dir)
        cfg = (run_dir / "solverConfig.yaml").read_text()
        note = run_dir / RA.MAIN_STAGE_STARTED
        seen.append({"cfg": cfg, "note": json.loads(note.read_text()) if note.exists() else None,
                     "init_sha": hashlib.sha256((run_dir / "nozzle.h5").read_bytes()).hexdigest()})
        n = RA._cfg_value(cfg, "nStepOuter")
        (run_dir / f"res_{n}.h5").write_bytes(f"field after {n}".encode())
        return 0

    def fake_restart(res_h5, mesh_h5):
        Path(mesh_h5).write_bytes(Path(res_h5).read_bytes())   # 最終場を初期の場へ (同じ格子の restart の代わり)

    RA.run_forge, RA.stage_gate, RA._restart_same_mesh = fake_forge, (lambda res, cfg: []), fake_restart
    try:
        for tag, cfg_m, cfg_p in (("ml", main, base), ("orig", main0, None)):
            seen.clear()
            rd = TMP / f"run_{tag}"
            rd.mkdir()
            (rd / "solverConfig.yaml").write_text(cfg_m)
            if cfg_p is not None:
                (rd / RA.PRESTAGE_CONFIG).write_text(cfg_p)
            (rd / "bcondConfig.yaml").write_text("wall: {physID: 3}\n")
            (rd / "nozzle.h5").write_bytes(b"initial field")
            (rd / RA.MAIN_STAGE_STARTED).write_text("{}")          # 前の呼び出しの通知 (最初に消えること)
            rc = RA.run_staged_ns(rd, stages="full")
            want_cfgs = [c for _, c, *_ in (pre if cfg_p is not None else pre0)] + [cfg_m]
            check(f"(d) {tag}: forge に渡った設定 = S1・S2・本段 (ns_stage_configs と文字列一致)",
                  rc == 0 and [s["cfg"] for s in seen] == want_cfgs)
            check(f"(d) {tag}: 前の段の間は通知が無く (古い通知は消えた)、本段の起動の時点で通知がある",
                  seen[0]["note"] is None and seen[1]["note"] is None and seen[2]["note"] is not None)
            note = seen[2]["note"]
            check(f"(d) {tag}: 通知の本段の設定の sha256 = 本段の solverConfig.yaml",
                  note["main_config"]["sha256"] == sha(cfg_m) and note["main_config"]["nStepOuter"] == 12000)
            check(f"(d) {tag}: 通知の初期の場の sha256 = 本段の起動の時点の nozzle.h5 (S2 の最終場を移した後)",
                  note["initial_field"] == {"file": "nozzle.h5", "sha256": seen[2]["init_sha"]}
                  and seen[2]["init_sha"] == hashlib.sha256(b"field after 3000").hexdigest())
            check(f"(d) {tag}: 通知の前の段の sha256 = forge に渡った S1・S2、前の段の出力は消えている",
                  [(s["tag"], s["nsteps"], s["config_sha256"]) for s in note["previous_stages"]]
                  == [("S1_soft", 3000, sha(seen[0]["cfg"])), ("S2_mid", 3000, sha(seen[1]["cfg"]))]
                  and note["res_files_present"] == [] and note["stages"] == "full" and "time" in note)
            check(f"(d) {tag}: 通知の前の段の元の設定 (main_line の run だけ)",
                  note["prestage_config"] == (None if cfg_p is None else {"file": RA.PRESTAGE_CONFIG, "sha256": sha(cfg_p)}))
            rec = RA.ns_stage_sha256(cfg_m, cfg_p)
            check(f"(d) {tag}: prepare_info の ns_stage_configs の sha256 = 実際に回した S1・S2・本段",
                  [rec["S1_soft"], rec["S2_mid"], rec["main"]] == [sha(s["cfg"]) for s in seen]
                  and rec["prestage_config"] == (None if cfg_p is None else sha(cfg_p)))
            man = json.loads((rd / "stage_manifest.json").read_text())
            check(f"(d) {tag}: stage_manifest の段は S1_soft・S2_mid・main", [s["tag"] for s in man["stages"]] == ["S1_soft", "S2_mid", "main"])
            if cfg_p is not None:
                check("(d) ml: 本段の区間は前の段と別の段の鍵 (check_convergence --segment が本段だけを取る)",
                      man["stages"][2]["key"] != man["stages"][1]["key"])
    finally:
        RA.run_forge, RA.stage_gate, RA._restart_same_mesh = orig

    # --- (e) 拒否 -----------------------------------------------------------------------------------------------------
    rejects("(e) main_line + evaluate.cfl_main", lambda: RA._main_line_spec(problem("e1", with_ml(ML, drop_cfl=False))), "cfl_main")
    rejects("(e) main_line + evaluate.implicit_relax の違う値", lambda: RA._main_line_spec(problem("e2", with_ml(ML, implicit_relax=0.5))),
            "implicit_relax")
    rejects("(e) main_line + 引数 cfl_main", lambda: RA._main_line_spec(p_ml, cfl_main=4.0), "引数 cfl_main")
    rejects("(e) main_line + 引数 implicit_relax の違う値", lambda: RA._main_line_spec(p_ml, implicit_relax=0.5), "implicit_relax")
    # 2026-10-13 親の決定: 前の段にも効く implicit_relax は main_line.relax と同じ値なら併記してよい (生産の系譜 run_0182 は S1・S2 も 0.7)
    for name, pr, kw in (("evaluate.implicit_relax 0.7", problem("e_r07", with_ml(ML, implicit_relax=0.7)), {}),
                         ("引数 implicit_relax 0.7", p_ml, {"implicit_relax": 0.7})):
        ml_r = RA._main_line_spec(pr, **kw)
        b_r, m_r = configs(pr, "e_r07_" + ("ev" if not kw else "arg"), ml=ml_r, **kw)
        check(f"(e) {name}: 受け付け、前の段の元の設定に implicitRelax 0.7", get(b_r, ("time", "deltaT", "implicitRelax")) == 0.7)
        check(f"(e) {name}: 本段の設定は relax を併記しないときと文字列一致", m_r == main, None if m_r == main else "差あり")
        pre_r, _, _ = RA.ns_stage_configs(m_r, "full", cfg_pre=b_r)
        check(f"(e) {name}: S1・S2 も implicitRelax 0.7 (B0 だけのキーは無い)",
              all(get(c, ("time", "deltaT", "implicitRelax")) == 0.7 and not [pth for pth in B0_ONLY if get(c, pth) is not RA._CFG_MISSING]
                  for _, c, _, _ in pre_r))
    for i, (name, mod) in enumerate((("line: only", {"line": "only"}), ("未知のキー isp", {"isp": 1}), ("itj 4 (ビット 1 無し)", {"itj": 4}),
                      ("itj true", {"itj": True}), ("itj 8", {"itj": 8}), ("relax 1.5", {"relax": 1.5}), ("cap 0", {"cap": 0}),
                      ("cfl 文字列", {"cfl": "4"}), ("limiter_ref に a が無い", {"limiter_ref": {"ro": 1.0, "p": 1.0}}),
                      ("limiter_ref の負の値", {"limiter_ref": {"ro": 1.0, "p": -1.0, "a": 1.0}}),
                      ("extra_fields の不正な名前", {"extra_fields": ["res_ro", "a, b"]}),
                      ("extra_fields の重複", {"extra_fields": ["res_ro", "res_ro"]}))):
        rejects(f"(e) {name}", lambda mod=mod, i=i: RA._main_line_spec(problem(f"e_bad{i}", with_ml({**ML, **mod}))))
    # limiter_ref の省略 = 自動 (plan tooling-settle-convergence-and-monitor §4.6、2026-10-16 に「欠けは拒否」から変更): 本段の設定に
    # 基準値を書かない (forge が本段の開始の場で決めて limiter_ref.json に書く)。それ以外は同じ
    src_auto = {}
    ml_auto = RA._main_line_spec(problem("e_nolr", with_ml({k: v for k, v in ML.items() if k != "limiter_ref"})), source_out=src_auto)
    m_auto = RA._apply_main_line(base, ml_auto)
    check("(e) limiter_ref の省略は自動 (limiter_ref {}・由来 auto; 戻りのキーは今までと同じ)", ml_auto["limiter_ref"] == {}
          and src_auto["limiter_ref"]["mode"] == "auto" and set(ml_auto) == set(ml))
    check("(e) 自動の本段の設定: 基準値のキーが無く、明示の本段との差は基準値 3 つだけ",
          all(get(m_auto, ("space", k)) is RA._CFG_MISSING for k in ("limiterRoRef", "limiterPRef", "limiterARef", "limiterRefLength"))
          and set(RA._cfg_diff_paths(RA._cfg_load(m_auto), RA._cfg_load(main)))
          == {("space", k) for k in ("limiterRoRef", "limiterPRef", "limiterARef")})
    rejects("(e) lowMachPrecond 2 の設定に B0", lambda: RA._apply_main_line(base.replace("lowMachPrecond: 0", "lowMachPrecond: 2"), ml),
            "lowMachPrecond")
    rejects("(e) 既にライン陰解法がある設定に B0", lambda: RA._apply_main_line(main, ml), "既に")
    rejects("(e) 前の段の元の設定がある run の stages ramp", lambda: RA.ns_stage_configs(main, "ramp", cfg_pre=base), "ramp")
    rejects("(e) 前の段の元の設定と本段の差が B0 以外 (convMethod)",
            lambda: RA.ns_stage_configs(main.replace("convMethod: 1", "convMethod: 2"), "full", cfg_pre=base), "convMethod")
    m2 = RA._apply_main_line(base, {**ml, "extra_fields": None})
    check("(e') extra_fields 無しなら output を書かない (他は同じ)",
          get(m2, ("output",)) is RA._CFG_MISSING and RA._cfg_diff_paths(RA._cfg_load(m2), RA._cfg_load(main)).keys() == {("output",)})
    pr = RA.ns_stage_configs(main.replace("nStepOuter: 12000", "nStepOuter: 200000"), "full", cfg_pre=base)
    check("(e') 本段の長さ (nStepOuter) を後で変えても前の段は同じ", [c for _, c, *_ in pr[0]] == [c for _, c, *_ in pre])

for s in SKIPPED:
    print("skip " + s)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILURES'}" + (f" ({len(SKIPPED)} skipped)" if SKIPPED else ""))
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
