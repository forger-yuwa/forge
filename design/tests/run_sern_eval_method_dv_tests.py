#!/usr/bin/env python3
"""評価方式 ID から設計変数を分離する試験 (plan tooling-sern-te-wake-grid §5.1 R3、codex result レビュー 2026-10-11 M3)。
forge は回さない (変換器があれば粗い 3D 格子を runner の prepare で作る)。

合格 (測る前に固定):
  (a) 方式 ID は設計変数の値に依らない: run_1084 (L_sw 0.8) と run_1096 (L_sw 1.0) の問題 YAML (同じ g3・局所化) の要求する
      評価方式の id が同じ (方式レコードが全項目で一致)。driver_sern の dv (M_c・f・θ_r0・θ_c0・L_cowl) と設計の角度でも変わらない
  (b) 方式に残るもの: L_sw の置き方 (L_sw_exact)・L_sw の指定の有無 (None = 側壁はカウル後縁まで)・解像度 (層数・第一層)・
      後縁下流の格子変形 (L_b 0 / band / local)・変換器の版で id が変わる
  (c) 設計変数の値は設計点・来歴に残る: prepare_info の eval_method.design_values (L_sw・設計の角度)・metrics の eval_method・
      台帳の作動点要約。粗い 3D 格子で L_sw 0.8 / 1.0 を prepare すると id が同じ (= キャンペーンの要求)・格子署名は違い、
      同じ系列のキャンペーンで両方の行が学習に入る。旧格子 (te_wake 0・band・旧変換器の sha・R3 前の L_sw 入りの id) は外れる

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_eval_method_dv_tests.py
"""
import json
import shutil
import sys
import tempfile
import types
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(TOOLS))
try:
    import pymoo  # noqa: F401
except ImportError:
    _m = types.ModuleType("forge_design.opt.moo"); _m.propose_infill = None
    sys.modules["forge_design.opt.moo"] = _m
try:
    import smt  # noqa: F401
except ImportError:
    _s = types.ModuleType("forge_design.opt.surrogate"); _s.KrigingSet = object
    sys.modules["forge_design.opt.surrogate"] = _s

import stage_manifest as SM  # noqa: E402
from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
from forge_design.probdef import Problem, load_problem  # noqa: E402

C46 = ROOT / "case" / "46.sern_design"
Y08 = C46 / "problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml"     # run_1084 の問題 YAML (L_sw 0.8、g3、局所化)
Y10 = C46 / "problem_3d_prod_3op_wallres_lswx10_tewake_g3_L10.yaml"     # run_1096 の問題 YAML (L_sw 1.0)
YA0 = C46 / "problem_3d_prod_3op_wallres_lswx08_tewake_g3_A0.yaml"      # 旧格子: te_wake 0
YB10 = C46 / "problem_3d_prod_3op_wallres_lswx08_tewake_g3_B10.yaml"    # 旧格子: band (変位の局所化の前)
PROD08, PROD10 = C46 / "problem_3d_prod_3op_wallres_lswx08.yaml", C46 / "problem_3d_prod_3op_wallres_lswx10.yaml"
fails = 0
skipped = []
tmp = Path(tempfile.mkdtemp(prefix="eval_method_dv_"))
CONV = {"sha256": "c" * 64}


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""))


def with_raw(p_, mesh=None, mesh3d=None, dv=None, drop3d=()):
    raw = json.loads(json.dumps(p_.raw))
    raw["mesh"].update(mesh or {}); raw.setdefault("mesh3d", {}).update(mesh3d or {}); raw["dv"].update(dv or {})
    for k in drop3d:
        raw["mesh3d"].pop(k, None)
    return Problem(name=p_.name, type=p_.type, gamma=p_.gamma, cp=p_.cp, spec=raw["spec"], dv=raw["dv"], geometry=raw["geometry"],
                   mesh=raw["mesh"], evaluate=raw["evaluate"], raw=raw)


def req(p_, conv=CONV):
    """キャンペーンの要求 (eval_method_required) を固定の変換器の版で。"""
    return R2.eval_method_required(p_, 3, converter=conv)


# ======================================================================== (a) 設計変数に依らない
p08, p10 = load_problem(Y08), load_problem(Y10)
r08, r10 = req(p08), req(p10)
check("(a) run_1084 (L_sw 0.8) と run_1096 (L_sw 1.0) の問題 YAML: 要求する評価方式の id が同じ",
      r08[R2.EVAL_METHOD_ID] == r10[R2.EVAL_METHOD_ID] is not None, (r08[R2.EVAL_METHOD_ID], r10[R2.EVAL_METHOD_ID]))
check("(a) 方式レコードが全項目で一致 (レビューで違っていた mesh_recipe_sha256 も同じ)・格子変形の実効値も同じ",
      r08["eval_method"] == r10["eval_method"] and r08[R2.TE_WAKE_EFFECTIVE] == r10[R2.TE_WAKE_EFFECTIVE] != "off")
check("(a) 生産 YAML (lswx08・lswx10) も同じ id (検証の YAML と name 以外は同一)",
      req(load_problem(PROD08))[R2.EVAL_METHOD_ID] == req(load_problem(PROD10))[R2.EVAL_METHOD_ID] == r08[R2.EVAL_METHOD_ID])
prm08, prm10 = R3.sern_mesh3d_params(p08, -0.70, 0.05), R3.sern_mesh3d_params(p10, -0.70, 0.05)
rc08, rc10 = R2.mesh_recipe(prm08, p08.mesh), R2.mesh_recipe(prm10, p10.mesh)
check("(a) レシピでは L_sw の値を 'design' に置き換え、他の項目は同じ", rc08["L_sw"] == rc10["L_sw"] == "design" and rc08 == rc10)
dv08, dv10 = R2.design_mesh_values(prm08, p08.mesh), R2.design_mesh_values(prm10, p10.mesh)
check("(a) 置き換えた値は design_mesh_values に残る (L_sw 0.8 / 1.0 と設計の角度)",
      dv08.get("L_sw") == 0.8 and dv10.get("L_sw") == 1.0 and set(dv08) == {"L_sw", "interface_angle", "top_ext_angle"}
      and dv08["interface_angle"] == -0.70, (dv08, dv10))
for k, v in (("M_c", 3.1), ("f", 0.5), ("theta_r0_deg", 17.0), ("theta_c0_deg", 4.0), ("L_cowl", 1.37)):
    check(f"(a) driver_sern の dv {k} を変えても id は同じ", req(with_raw(p08, dv={k: v}))[R2.EVAL_METHOD_ID] == r08[R2.EVAL_METHOD_ID])
eid = lambda p_, thb, the: R2.eval_method_of(3, "node", R2.mesh_recipe(R3.sern_mesh3d_params(p_, thb, the), p_.mesh), R2.mesher_version(3),
                                            *R2.te_wake_spec(p_.raw, 3), CONV["sha256"], R2.te_wake_mode_spec(p_.raw, 3))[1]
check("(a) 設計の角度 (θ_b・θ_e) を変えても id は同じ、要求 (仮の角度) とも同じ",
      eid(p08, -0.70, 0.05) == eid(p08, -0.50, 0.10) == eid(p10, -0.60, 0.0) == r08[R2.EVAL_METHOD_ID])
check("(a) L_sw の値を他の値 (0.9) にしても id は同じ", req(with_raw(p08, mesh3d={"L_sw": 0.9}))[R2.EVAL_METHOD_ID] == r08[R2.EVAL_METHOD_ID])

# ======================================================================== (b) 方式に残るもの
for label, kw in (("L_sw_exact true → false", {"mesh3d": {"L_sw_exact": False}}),
                  ("L_sw の指定なし (側壁はカウル後縁まで)", {"drop3d": ("L_sw",)}),
                  ("層数 nz_in", {"mesh3d": {"nz_in": 26}}), ("層数 nj_bot", {"mesh3d": {"nj_bot": 56}}),
                  ("第一層 first_wall_frac", {"mesh3d": {"first_wall_frac": 2.0e-4}}), ("ni_noz", {"mesh3d": {"ni_noz": 91}}),
                  ("局所化 → band (te_wake_mode)", {"drop3d": ("te_wake_mode",)}), ("te_wake 0 (旧格子)", {"mesh3d": {"te_wake_blend_H": 0.0}}),
                  ("L_b 0.5", {"mesh3d": {"te_wake_blend_H": 0.5}})):
    check(f"(b) {label} で id が変わる", req(with_raw(p08, **kw))[R2.EVAL_METHOD_ID] != r08[R2.EVAL_METHOD_ID])
check("(b) 旧格子の問題 YAML (te_wake 0 = A0・band = B10) は別の id", len({r08[R2.EVAL_METHOD_ID], req(load_problem(YA0))[R2.EVAL_METHOD_ID],
                                                                   req(load_problem(YB10))[R2.EVAL_METHOD_ID]}) == 3)
check("(b) 変換器の版違いで id が変わる", req(p08, {"sha256": "d" * 64})[R2.EVAL_METHOD_ID] != r08[R2.EVAL_METHOD_ID])
check("(b) L_sw の指定なしのレシピは None のまま (方式として残る)、design_values に L_sw は無い",
      (lambda p_: R2.mesh_recipe(R3.sern_mesh3d_params(p_, 0.0, 0.0), p_.mesh)["L_sw"] is None
       and "L_sw" not in R2.design_mesh_values(R3.sern_mesh3d_params(p_, 0.0, 0.0), p_.mesh))(with_raw(p08, drop3d=("L_sw",))))
# R3 の前の定義 (L_sw の値がレシピに入る) で作った id は、新しい要求と一致しない (旧 run の記録は再判定で外れる)
_old_recipe = {**rc08, "L_sw": 0.8}
_old_id = R2.eval_method_of(3, "node", _old_recipe, R2.mesher_version(3), *R2.te_wake_spec(p08.raw, 3), CONV["sha256"],
                            R2.te_wake_mode_spec(p08.raw, 3))[1]
check("(b) R3 の前の定義 (L_sw の値入りのレシピ) の id は新しい要求と別", _old_id != r08[R2.EVAL_METHOD_ID])

# ======================================================================== (c) 実際の prepare (変換器があるとき)
HAVE_CONV = R2.CONVERTER.exists()


class Fake:
    def __init__(self, rows, te="off", req_id=None, co=()):
        self.rows, self.te_wake_required, self.eval_method_required, self.eval_method_carry_over = rows, te, req_id, tuple(co)
        self.ref, self.ops = (-0.9, 20.0), []


def stub_collect(rd, prob_path):
    """力の履歴とゲートを差し替えて 3D の collect (metrics.json) を回す。起動記録は本段の config で 1 行。"""
    cfg_main = (rd / "solverConfig_main.yaml").read_text()
    with open(rd / "forge_launches.jsonl", "a") as fh:
        fh.write(json.dumps({"time": 0, "cfg_fnv": SM.fnv1a64(cfg_main), "slauWallNormalChi": 1, "scalarGradient": "lsq"}) + "\n")
    gates = lambda *a, **k: {"steadiness": {"series": {}}, "objective": "C_T", "verdict": "PASS", "fail_class": None, "reasons": []}
    saved = R3.evaluate_gates
    R3.evaluate_gates = gates
    try:
        return R3.collect(prob_path, rd, rc=0)
    finally:
        R3.evaluate_gates = saved


COARSE3D = {"ni_up": 6, "ni_noz": 30, "ni_plume": 40, "nj_top": 21, "nj_bot": 15, "nz_in": 4, "nz_out": 3, "nj_vside": 5,
            "first_wall_frac": 2.0e-3, "first_z_frac": 0.02, "first_wake_frac": 0.004}
COARSE2 = {"nj_ext_top": 9, "first_top_frac": 0.004}
if not HAVE_CONV:
    skipped.append("変換器が無い: 実際の prepare の試験")
    print(f"SKIP (c) 実際の prepare: {R2.CONVERTER} が無い")
else:
    cache = {}
    _saved_design = R2.design_from_problem

    def _cached_design(p_, design=None):
        key = json.dumps(p_.dv, sort_keys=True, default=str)
        if key not in cache:
            cache[key] = _saved_design(p_, design=design)
        return cache[key]
    R2.design_from_problem = _cached_design
    runs = {}
    try:
        for lab, src in (("L08", Y08), ("L10", Y10), ("A0", YA0), ("B10", YB10)):
            y = yaml.safe_load(open(src))
            y["mesh3d"].update(COARSE3D); y["mesh"].update(COARSE2)
            pth = tmp / f"p_{lab}.yaml"; pth.write_text(yaml.safe_dump(y, allow_unicode=True))
            runs[lab] = (tmp / f"r_{lab}", pth, R3.prepare(pth, tmp / f"r_{lab}", op="m10_on"))
    finally:
        R2.design_from_problem = _saved_design
    (d08, q08, i08), (d10, q10, i10) = runs["L08"], runs["L10"]
    e08, e10 = i08["eval_method"], i10["eval_method"]
    rq = R2.eval_method_required(load_problem(q08), 3)
    check("(c) 粗い 3D の L_sw 0.8 / 1.0: prepare の id が同じ = キャンペーンの要求 (問題 YAML・現在のコード・現在の変換器)",
          e08[R2.EVAL_METHOD_ID] == e10[R2.EVAL_METHOD_ID] == rq[R2.EVAL_METHOD_ID] is not None)
    check("(c) 格子署名・双対幾何のハッシュ・実際の側壁の後縁 (minfo の L_sw) は違う",
          i08["mesh_provenance"]["mesh_signature"] != i10["mesh_provenance"]["mesh_signature"]
          and i08["mesh_provenance"]["dual_hash"] != i10["mesh_provenance"]["dual_hash"] and i08["mesh"]["L_sw"] != i10["mesh"]["L_sw"])
    check("(c) prepare_info の eval_method.design_values に L_sw 0.8 / 1.0 と設計の角度 (prepare の θ_b と一致)",
          e08["design_values"]["L_sw"] == 0.8 and e10["design_values"]["L_sw"] == 1.0
          and abs(e08["design_values"]["interface_angle"] - __import__("math").radians(i08["design"]["theta_b_deg"])) < 1e-12
          and e08["mesh_recipe"]["L_sw"] == "design", e08.get("design_values"))
    outs = {lab: stub_collect(runs[lab][0], runs[lab][1]) for lab in runs}
    m08 = json.loads((d08 / "metrics.json").read_text())
    check("(c) metrics.json の eval_method に design_values (L_sw 0.8) と id",
          m08["eval_method"]["design_values"]["L_sw"] == 0.8 and m08[R2.EVAL_METHOD_ID] == e08[R2.EVAL_METHOD_ID])
    summ = {lab: {**D.SernCampaign._op_summary(outs[lab], runs[lab][0]), "scalar_gradient_effective": "lsq", "slau_wall_normal_chi_effective": 1}
            for lab in runs}
    check("(c) 台帳の作動点要約に design_values (設計点) が残る",
          summ["L08"]["design_values"]["L_sw"] == 0.8 and summ["L10"]["design_values"]["L_sw"] == 1.0)

    def row_(s, x):
        return {"status": "PASS", "flag_policy": D.FLAG_POLICY, "x": [x], "C_T_w": 0.9, "L_ramp": 1.0, "C_M_w": 0.0, "ops": {"m10_on": s}}
    TE = rq[R2.TE_WAKE_EFFECTIVE]
    learn = lambda rows: len(D.SernCampaign._XF(Fake(rows, TE, rq[R2.EVAL_METHOD_ID]))[0])
    check("(c) L_sw 0.8 / 1.0 の 2 行は同じ系列 (局所化のキャンペーン) で一緒に学習に入る",
          learn([row_(summ["L08"], 0.8), row_(summ["L10"], 1.0)]) == 2)
    check("(c) 旧格子 te_wake 0 (A0)・band (B10) の行は外す (id が違う)",
          learn([row_(summ["A0"], 0.0)]) == 0 and learn([row_(summ["B10"], 0.0)]) == 0
          and len({summ[k][D.EVAL_METHOD_ID] for k in ("L08", "A0", "B10")}) == 3)
    check("(c) 旧変換器の sha の行・R3 の前の定義の id (L_sw 入り) の行は外す",
          learn([row_({**summ["L08"], D.EVAL_METHOD_ID: R2.eval_method_of(3, "node", e08["mesh_recipe"], e08["mesher"], 1.0,
                                                                         e08["eval_method"]["te_wake_curve_version"], "d" * 64,
                                                                         e08["eval_method"].get("te_wake_mode_version"))[1]}, 0.8)]) == 0
          and learn([row_({**summ["L08"], D.EVAL_METHOD_ID: R2.eval_method_of(
              3, "node", {**e08["mesh_recipe"], "L_sw": 0.8}, e08["mesher"], 1.0, e08["eval_method"]["te_wake_curve_version"],
              e08["converter"]["sha256"], e08["eval_method"].get("te_wake_mode_version"))[1]}, 0.8)]) == 0)

shutil.rmtree(tmp, ignore_errors=True)
if skipped:
    print("省略: " + "; ".join(skipped))
print("VERDICT:", ("PASS" if not skipped else "PASS (一部省略)") if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
