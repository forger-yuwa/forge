"""成果物の一括書き出し `forge_design.export.deliverables` の試験 (plan tooling-nozzle-wt-deliverables §6 V6・V7 のゲート部分)。

fixture は FORGE_DELIV_FIXTURE (既定: このセッションの scratchpad の m6_prod; 作り方は tests/make_deliverables_fixture.py)。
M6 の生産の壁の準備ディレクトリ (wall_repr.json・prepare_info.json・bcondConfig.yaml・nozzle.msh) を仮の run にして、
(a) 一括書き出しが通り MANIFEST の sha256 が実ファイルと一致、(b) 寸法ゲート: 生産の設計値で作ったカードは PASS・出口長を 2.2 mm
ずらしたカードは FAIL で MANIFEST の verdict も FAIL、(c) 入力が欠けた run (wall_repr.json なし) は cad を理由つきで飛ばす、を確かめる。
FreeCAD・gmsh が要る (無ければ飛ばす)。作業は一時ディレクトリで行い、fixture には書かない。

壁法線の近壁層の格子 (構造ファイル v2、plan tooling-nozzle-card-local-mesh §5.1 #5, 2026-10-15):
- (S) fixture なしで回る: mesh2d の小さな格子の run で part_mesh が、縦線の格子は v1 (変更前 = git の BASE_COMMIT の
  mesh_structure.py の出力とバイト一致)、壁法線の近壁層の格子で壁ファイルが無い run は「作らない」と理由 (壁の傾き) を書く。
- (W) FORGE_DELIV_FIXTURE_WN (壁法線の近壁層の run から make_deliverables_fixture.py で作ったもの。既定は下の FIXTURE_WN_DEFAULT =
  case/71 run_0003_ns_cont の写し) があれば、part_mesh が v2 を作り STEP からの作り直しの比較が PASS。格子品質は変換器が
  FP64 の座標でないと薄い第 1 層を測れない (FORGE_CONVERTER で指定) ので、元と作り直しの判定が同じことだけを合否にする。
  規則の照合 (prepare_info.json の mesh.params の wall_normal_layer・wall_normal_end_stations): part_mesh で実施されて通り、
  同じ格子で傾きを (1 + 1e-6) 倍・d_n を 1 % ずらすと StructureError、パラメータ無しでは通って「未実施」と書く。
"""
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.export import deliverables as dl  # noqa: E402

FIX = Path(os.environ.get("FORGE_DELIV_FIXTURE",
                          "/tmp/claude-1000/-home-sano-work-forge/bd18adac-85f4-4627-811b-8f3186867995/scratchpad/fixtures/m6_prod"))
CFG_TEMPLATE = Path(os.environ.get("FORGE_DELIV_CFG_TEMPLATE",
                                   "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0051_ns_final_c2/solverConfig.yaml"))
FIXTURE_WN_DEFAULT = "/tmp/claude-1000/-home-sano-work-forge/e254baf2-8857-4243-99d2-53804f553ad9/scratchpad/m55_deliv/fixtures/m55_wn"
FIX_WN = Path(os.environ.get("FORGE_DELIV_FIXTURE_WN", FIXTURE_WN_DEFAULT))
BASE_COMMIT = "6d4334e3"     # v2 を入れる前の mesh_structure.py (v1 の出力のビット同一の基準)
fails = 0


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info else "")); fails += (not ok)


def make_run(tmp: Path, with_wall=True) -> Path:
    run = tmp / "run_fixture"
    run.mkdir()
    for n in ("prepare_info.json", "bcondConfig.yaml", "nozzle.msh") + (("wall_repr.json",) if with_wall else ()):
        os.symlink(FIX / n, run / n)       # 大きいファイルはリンク (deliverables は run に書かない; 出力は --out)
    if CFG_TEMPLATE.exists():
        shutil.copy2(CFG_TEMPLATE, run / "solverConfig.yaml")
    return run


def card_from_production(L_div_mm: float, throat_d_mm: float) -> dict:
    rt = 0.07666536551630307
    return {"format": "forge.wt_card/v1", "inlet": {"radius_m": 6.485 * rt, "pipe_length_m": 0.5 * rt, "contraction_length_m": 12.0 * rt},
            "divergent_length_m": L_div_mm / 1000.0, "throat": {"diameter_physical_m": throat_d_mm / 1000.0},
            "tolerances": {"length_m": 1.0e-3, "diameter_m": 2.0e-5}}


def base_mesh_structure(tmp: Path):
    """BASE_COMMIT の mesh_structure.py を git から読んで import する (読めなければ None)。"""
    import importlib.util
    import subprocess
    r = subprocess.run(["git", "-C", str(Path(__file__).resolve().parents[2]), "show",
                        f"{BASE_COMMIT}:design/forge_design/export/mesh_structure.py"], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    p = tmp / "mesh_structure_base.py"
    p.write_text(r.stdout, encoding="utf-8")
    spec = importlib.util.spec_from_file_location("mesh_structure_base", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def structure_tests() -> None:
    """(S) part_mesh の構造ファイル: 縦線の格子は v1 (変更前とバイト一致)、壁ファイルの無い壁法線の近壁層の格子は理由つきで作らない。"""
    import numpy as np
    from forge_design.meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
    from forge_design.export import mesh_structure as ms
    if not Path(dl.GMSH_PYTHON).exists():
        print(f"skip (S): gmsh の python が無い ({dl.GMSH_PYTHON})")
        return

    class WallS:
        x_in, x_e = -2.0, 4.0

        @staticmethod
        def r(x, deriv=0):
            x = np.asarray(x, dtype=float)
            return 1.0 + 0.05 * x * x if deriv == 0 else 0.1 * x

    with tempfile.TemporaryDirectory(prefix="deliv_struct_") as td:
        tmp = Path(td)
        base = base_mesh_structure(tmp)
        for name, extra in (("v", {}), ("w", {"wall_normal_layer": [0.02, 0.8]})):
            run = tmp / f"run_{name}"
            run.mkdir()
            c, q, b = generate_axisym_mesh(WallS, Mesh2DParams(ni=41, nj=17, wall_first_frac=0.01, scale=0.05, **extra))
            write_msh41_2d(run / "nozzle.msh", c, q, b, digits=17)
            res = dl.part_mesh(run, tmp / f"out_{name}")
            if name == "v":
                js = tmp / "out_v" / "mesh" / "mesh_structure.json"
                check("(S1) 縦線の格子: part_mesh が v1 の構造ファイルを書く (STEP が無いので作り直しは飛ばす)",
                      res.get("structure_format") == ms.FORMAT and js.exists() and "STEP が無い" in str(res.get("rebuild")), res.get("structure"))
                if base is None:
                    print(f"skip (S1) v1 のビット同一: git から {BASE_COMMIT} を読めない")
                else:
                    cc, qq, bb, _ = base.read_msh41_mesh2d(tmp / "out_v" / "mesh" / "nozzle.msh")
                    base.write_structure_json(base.extract_structure(cc, qq, bb), tmp / "base.json")
                    check(f"(S1) 縦線の格子の mesh_structure.json が変更前 ({BASE_COMMIT} の extract_structure) とバイト一致",
                          js.exists() and js.read_bytes() == (tmp / "base.json").read_bytes())
            else:
                st = str(res.get("structure"))
                check("(S2) 壁法線の近壁層の格子で壁ファイルの無い run: 構造ファイルを作らず理由 (壁の傾きを渡す・壁の傾きを作れない) を書く",
                      st.startswith("作らない") and "壁の傾きを渡す" in st and "壁の傾きを作れない" in st
                      and not (tmp / "out_w" / "mesh" / "mesh_structure.json").exists(), st[-160:])


def wn_fixture_tests() -> None:
    """(W) 壁法線の近壁層の run の fixture (FORGE_DELIV_FIXTURE_WN) で part_mesh が v2 を作り、STEP からの作り直しの比較が PASS。"""
    import numpy as np
    from forge_design.export import mesh_structure as ms
    need = ("nozzle.msh", "wall_repr.json", "prepare_info.json", "bcondConfig.yaml")
    if not all((FIX_WN / n).exists() for n in need) or not (FIX_WN / "step" / "wall_physical.step").exists():
        print(f"skip (W): 壁法線の近壁層の fixture が無い ({FIX_WN})。tests/make_deliverables_fixture.py で作り FORGE_DELIV_FIXTURE_WN で指定")
        return
    with tempfile.TemporaryDirectory(prefix="deliv_wn_") as td:
        tmp = Path(td)
        run = tmp / "run_wn"
        run.mkdir()
        for n in need + ("solverConfig.yaml",):
            if (FIX_WN / n).exists():
                os.symlink(FIX_WN / n, run / n)
        out = tmp / "out"
        (out / "cad").mkdir(parents=True)
        shutil.copy2(FIX_WN / "step" / "wall_physical.step", out / "cad" / "wall_physical.step")
        res = dl.part_mesh(run, out)
        sc = res.get("structure_self_check") or {}
        src = res.get("structure_wall_slope_source") or {}
        check("(W1) part_mesh が v2 (節点の規則 wall_normal_blend_v1) の構造ファイルを書く・自己検査 ≤ 1e-7 mm かつ第 1 層厚の 1e-3",
              res.get("structure_format") == ms.FORMAT_V2 and res.get("structure_node_rule") == ms.NODE_RULE_V2
              and sc.get("max_dist_mm", 1.0) <= 1e-7 and sc.get("max_dist_over_first_layer", 1.0) <= 1e-3,
              f"{res.get('structure')}: 自己検査 {sc.get('max_dist_mm', float('nan')):.2e} mm・{sc.get('max_dist_over_first_layer', float('nan')):.2e}")
        check("(W1) 壁の傾きの出典 = run の壁ファイル (sha256 一致)",
              src.get("file") == "wall_repr.json" and src.get("sha256") == hashlib.sha256((FIX_WN / "wall_repr.json").read_bytes()).hexdigest())
        rc = res.get("structure_rule_check") or {}
        rpar = res.get("structure_rule_params") or {}
        pinfo = json.loads((FIX_WN / "prepare_info.json").read_text())["mesh"]["params"]
        check("(W1) 規則の照合を実施して ≤ 1e-7 mm かつ第 1 層厚の 1e-3・パラメータは prepare_info.json の mesh.params (sha256・msh の sha256 一致)",
              rc.get("status") == "実施" and rc.get("max_dist_mm", 1.0) <= 1e-7 and rc.get("max_dist_over_first_layer", 1.0) <= 1e-3
              and rpar.get("wall_normal_layer") == [float(v) for v in pinfo["wall_normal_layer"]]
              and rpar.get("wall_normal_end_stations") == pinfo.get("wall_normal_end_stations")
              and rpar.get("source", {}).get("key") == "mesh.params" and rpar.get("source", {}).get("msh_sha256_matches") is True
              and rpar.get("source", {}).get("sha256") == hashlib.sha256((FIX_WN / "prepare_info.json").read_bytes()).hexdigest(),
              f"{rc.get('status')}: {rc.get('max_dist_mm', float('nan')):.2e} mm・{rc.get('max_dist_over_first_layer', float('nan')):.2e}・"
              f"|β_規則 − β_解| {rc.get('beta_rule_minus_solved_max_abs', float('nan')):.1e}・パラメータ {rpar.get('wall_normal_layer')}・"
              f"{rpar.get('wall_normal_end_stations')}")
        js = out / "mesh" / "mesh_structure.json"
        if js.exists():
            st = ms.load_structure_json(js)
            cc = ms.read_msh41_mesh2d(FIX_WN / "nozzle.msh")[0]
            e = float(np.abs(ms.structure_coords(st) - cc).max()) / (st["length_mm"] * 1e-3)
            check("(W2) 構造ファイルから作り直した座標が元の msh と ≤ 1e-12 L", e <= 1e-12, f"{e:.2e} L")
        rb = res.get("rebuild") if isinstance(res.get("rebuild"), dict) else {}
        bd = rb.get("b_distance", {})
        ec = rb.get("end_columns", {})
        check("(W3) STEP からの作り直しの比較が PASS (v2)", rb.get("verdict") == "PASS" and rb.get("structure_format") == ms.FORMAT_V2,
              f"最大距離 {bd.get('max_dist_mm', float('nan')):.3e} mm・第 1 層厚比 {bd.get('max_dist_over_first_layer', float('nan')):.3e}・"
              f"端の列の外れ 入口 {ec.get('inlet', {}).get('max_abs_dx_mm', float('nan')):.3e} mm・出口 "
              f"{ec.get('outlet', {}).get('max_abs_dx_mm', float('nan')):.3e} mm・理由 {rb.get('fail_reasons', str(res.get('rebuild'))[-200:])}")
        qo, qr = (res.get("quality", {}).get(k, {}) for k in ("original", "rebuilt"))
        check("(W4) 格子品質の判定が元と作り直しで同じ", qo.get("verdict") is not None and qo.get("verdict") == qr.get("verdict"),
              f"元 {qo.get('verdict')}・作り直し {qr.get('verdict')} (変換器 {qo.get('converter')}; FP64 でない変換器は薄い第 1 層で FAIL)")
        print(f"     (記録) part_mesh の pass = {res.get('pass')} (格子品質 元 {qo.get('verdict')}・作り直し {qr.get('verdict')})")
        # 実データでの規則の照合の負例 (同じ格子・同じ傾きの作り方で、傾きを定数倍・パラメータをずらす)
        cc, qq, bb, ph = ms.read_msh41_mesh2d(FIX_WN / "nozzle.msh")
        slope, _, why = dl._wall_slope(run, cc, bb)
        rule, why_r = dl._rule_params(run, hashlib.sha256((FIX_WN / "nozzle.msh").read_bytes()).hexdigest())
        if slope is None or rule is None:
            check("(W5) 実データの規則の照合の負例: 傾き・パラメータを作れる", False, f"{why} / {why_r}")
            return
        for name, kw, word in (("傾きを (1 + 1e-6) 倍", {"wall_slope": slope * (1.0 + 1e-6), "rule_params": rule}, "規則の照合"),
                               ("d_n を 1 % ずらす", {"wall_slope": slope, "rule_params": dict(
                                   rule, wall_normal_layer=[rule["wall_normal_layer"][0] * 1.01, rule["wall_normal_layer"][1]])}, "規則の照合")):
            try:
                ms.extract_structure(cc, qq, bb, phys=ph, **kw)
                check(f"(W5) 実データ: {name} と規則の照合で StructureError", False, "例外なし")
            except ms.StructureError as e:
                check(f"(W5) 実データ: {name} と規則の照合で StructureError", word in str(e), str(e)[:170])
        st0 = ms.extract_structure(cc, qq, bb, phys=ph, wall_slope=slope)
        check("(W5) 実データ: 規則のパラメータ無しでは今までどおり通り、checks に「未実施」",
              st0["checks"]["rule_check"]["status"] == "未実施 (パラメータ無し)" and st0["rule_params"] is None)
        st1 = ms.extract_structure(cc, qq, bb, phys=ph, wall_slope=slope * (1.0 + 1e-6))
        print(f"     (記録) 実データ: 傾きを (1 + 1e-6) 倍でも規則のパラメータ無しでは通る (自己検査 {st1['checks']['self_check']['max_dist_mm']:.2e} mm、"
              f"β が [0, 1] の外に出た節点 {st1['checks']['n_beta_clamped']})")


def main() -> int:
    structure_tests()
    wn_fixture_tests()
    if not (FIX / "wall_repr.json").exists() or not (FIX / "nozzle.msh").exists():
        print(f"skip: fixture が無い ({FIX})。tests/make_deliverables_fixture.py で作る")
        print(f"FAIL 件数: {fails}")
        return 1 if fails else 0
    with tempfile.TemporaryDirectory(prefix="deliv_test_") as td:
        tmp = Path(td)
        run = make_run(tmp)
        # 寸法は STEP から測る: 生産の設計値 (入口・直管・縮流部) と、STEP の終点・物理スロート径をカードにする
        from forge_design.export import dimension_gate as dg
        sd = dg.step_dimensions(FIX / "step" / "wall_physical.step")
        card_ok = card_from_production(sd["L_div_mm"], sd["throat_physical"]["diameter_mm"])
        cp = tmp / "card_ok.yaml"
        cp.write_text(yaml.safe_dump(card_ok))
        out = tmp / "deliv_ok"
        man = dl.export(run, out, card=cp)
        check("(a) 一括書き出しの verdict が PASS", man["verdict"] == "PASS",
              f"verdict {man['verdict']} failed {man['failed_parts']} errors {list(man['errors'])} skipped {list(man['skipped'])}")
        bad = [k for k, h in man["files_sha256"].items() if hashlib.sha256((out / k).read_bytes()).hexdigest() != h]
        check("(a) MANIFEST の sha256 が実ファイルと一致", not bad and len(man["files_sha256"]) > 5, f"{len(man['files_sha256'])} files, bad {bad}")
        for k in ("cad/wall_physical.step", "cad/wall_surface.step", "mesh/nozzle.msh", "mesh/mesh_structure.json",
                  "mesh/rebuild_from_step.py", "fluent/nozzle_2d_axi.msh", "fluent/boundary_conditions.md",
                  "excel/wall_shape_mm.csv", "DIMENSIONS.json", "README.md"):
            check(f"(a) {k} がある", (out / k).exists())
        base = base_mesh_structure(tmp)
        if base is not None and (out / "mesh" / "mesh_structure.json").exists():
            cc, qq, bb, _ = base.read_msh41_mesh2d(FIX / "nozzle.msh")
            base.write_structure_json(base.extract_structure(cc, qq, bb), tmp / "base_m6.json")
            check(f"(a) 縦線の格子 (M6) の mesh_structure.json が変更前 ({BASE_COMMIT}) とバイト一致",
                  (out / "mesh" / "mesh_structure.json").read_bytes() == (tmp / "base_m6.json").read_bytes())
        d = json.loads((out / "DIMENSIONS.json").read_text())
        check("(b) 寸法ゲート PASS (生産の設計値のカード)", d["verdict"] == "PASS",
              {k: round(v["diff_mm"], 6) for k, v in d.get("items", {}).items()})
        card_bad = dict(card_ok, divergent_length_m=card_ok["divergent_length_m"] + 0.0022)
        cb = tmp / "card_bad.yaml"
        cb.write_text(yaml.safe_dump(card_bad))
        out_b = tmp / "deliv_bad"
        man_b = dl.export(run, out_b, card=cb, skip=("fluent", "excel"))
        db = json.loads((out_b / "DIMENSIONS.json").read_text())
        check("(b) 出口長を 2.2 mm ずらしたカードは寸法ゲート FAIL・MANIFEST も FAIL",
              db["verdict"] == "FAIL" and man_b["verdict"] == "FAIL" and not db["items"]["L_div"]["pass"],
              f"L_div diff {db['items']['L_div']['diff_mm']:.4f} mm")
        check("(b) README の冒頭に寸法ゲートの FAIL が書かれる", "**FAIL**" in (out_b / "README.md").read_text().split("## 単位")[0])
        run2 = tmp / "nowall"
        run2.mkdir()
        for n in ("prepare_info.json", "bcondConfig.yaml", "nozzle.msh"):
            os.symlink(FIX / n, run2 / n)
        man_c = dl.export(run2, tmp / "deliv_nowall", card=cp, skip=("fluent",))
        check("(c) wall_repr.json の無い run は cad を理由つきで飛ばす", "cad" in man_c["skipped"] and "wall_repr" in man_c["skipped"]["cad"],
              man_c["skipped"].get("cad"))
        check("(c) その run の寸法ゲートは STEP が無いので飛ばす (合格にしない)", man_c["verdict"] != "PASS", man_c["verdict"])
    # ---- codex result レビュー 2026-10-13 の負例 (M1・M2・M4) ----
    import copy
    from unittest import mock
    from forge_design.export import wall_step
    md = dg.mesh_dimensions(FIX / "nozzle.msh")
    chm = dg.chord_from_mesh(FIX / "step" / "wall_physical.step", FIX / "nozzle.msh")   # 弦の検査は必須 (欠けたら未判定)
    check("(M1) 基準: 生産の設計値のカード・メッシュ・弦あり → PASS", dg.evaluate_gate(card_ok, sd, md, None, chm)["verdict"] == "PASS")
    for name, mut in (("入口半径 +10 mm", lambda m: m.__setitem__("inlet_radius_mm", m["inlet_radius_mm"] + 10.0)),
                      ("スロート半径 +10 mm", lambda m: m["throat_physical"].__setitem__("r_mm", m["throat_physical"]["r_mm"] + 10.0)),
                      ("直管の終端の区間を別の位置に", lambda m: m.__setitem__("straight_end_interval_mm", [0.0, 1.0])),
                      ("出口端 +0.01 mm", lambda m: m.__setitem__("x_end_mm", m["x_end_mm"] + 0.01))):
        m2 = copy.deepcopy(md); mut(m2)
        check(f"(M1) メッシュ側の {name} → 寸法ゲート FAIL", dg.evaluate_gate(card_ok, sd, m2, None, chm)["verdict"] == "FAIL")
    check("(M1) 弦の形状誤差 50 µm → FAIL",
          dg.evaluate_gate(card_ok, sd, md, None, {"status": "ok", "min_um": -50.0, "max_um": 3.0})["verdict"] == "FAIL")
    card_nt = {k: v for k, v in card_ok.items() if k != "tolerances"}
    g_nt = dg.evaluate_gate(card_nt, sd, md, None, chm)
    check("(M4) 公差を省いたカードは生成器と同じ既定の公差で判定 (止まらない)",
          g_nt["verdict"] == "PASS" and all(v.startswith("default") for v in g_nt["tolerances_source"].values()), g_nt["tolerances_source"])
    with tempfile.TemporaryDirectory(prefix="deliv_m2_") as td2:
        tmp2 = Path(td2); run3 = make_run(tmp2)
        real_ws = wall_step.write_surface
        def bad_surface(*a, **k):
            r = real_ws(*a, **k); r.update(readback_valid=False, readback_faces=2); return r
        with mock.patch.object(wall_step, "write_surface", bad_surface):
            pc = dl.part_cad(run3, tmp2 / "out")
        check("(M2) 3D STEP の読み直しが無効・面 2 枚なら cad は不合格", pc["pass"] is False and not pc["checks"]["freecad_readback_valid"]
              and not pc["checks"]["freecad_one_face"], {k: v for k, v in pc["checks"].items() if not v})
        with mock.patch.object(dl, "GMSH_PYTHON", "/bin/false"):
            pc2 = dl.part_cad(run3, tmp2 / "out2")
        check("(M2) gmsh の読み直しが失敗したら cad は不合格 (既定値で合格にしない)", pc2["pass"] is False and not pc2["checks"]["gmsh_readback"])
    # ---- codex result レビュー 2 回目の負例 ----
    md0 = dg.mesh_dimensions(FIX / "nozzle.msh")
    m3 = copy.deepcopy(md0); m3["x_start_mm"] += 10.0
    check("(R2-M1) メッシュの入口端を 10 mm ずらす (直管が短い) → 寸法ゲート FAIL", dg.evaluate_gate(card_ok, sd, m3, None, chm)["verdict"] == "FAIL")
    ch = dg.chord_from_mesh(FIX / "step" / "wall_physical.step", FIX / "nozzle.msh")
    check(f"(R2-M1) 弦の形状誤差を配布メッシュと STEP から測る (最大 {max(abs(ch['min_mm']), abs(ch['max_mm'])):.2e} mm) → 基準は PASS",
          dg.evaluate_gate(card_ok, sd, md0, None, ch)["verdict"] == "PASS")
    check("(R2-M1) 弦の検査が欠けたら PASS にしない (未判定)", dg.evaluate_gate(card_ok, sd, md0, None, None)["verdict"] == "未判定")
    check("(R2-M1) メッシュの測定が欠けたら PASS にしない (未判定)", dg.evaluate_gate(card_ok, sd, None, None, ch)["verdict"] == "未判定")
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_wall_tables_tests as RWT
    from forge_design.export import wall_tables as WT
    with tempfile.TemporaryDirectory(prefix="deliv_int_") as td3:
        syn = Path(td3) / "syn"
        RWT.make_synth(syn, steps=tuple(range(100, 2100, 100)), tw=1e-5, residual=True)
        okx = dl.part_excel(syn, Path(td3) / "o_ok")
        check(f"(R2-M2) 基準: 合成 run の積分の整合 PASS (相対 {okx['cfd_integral']['rel_to_sum_abs_Q']:.1e})・Excel part 合格",
              okx["cfd_integral"]["pass"] and okx.get("pass") is True)
        real_t = WT.cfd_wall_table
        def bad_table(run, out_csv, *a, **k):
            r = real_t(run, out_csv, *a, **k)
            lines = Path(out_csv).read_text().splitlines()
            hdr = lines[0].split(","); ic = hdr.index("ring_area_m2")
            out = [lines[0]] + [",".join(f"{float(v) * 1.01!r}" if j == ic else v for j, v in enumerate(ln.split(","))) for ln in lines[1:]]
            Path(out_csv).write_text("\n".join(out) + "\n")
            return r
        with mock.patch.object(WT, "cfd_wall_table", bad_table):
            man = dl.export(syn, Path(td3) / "o_bad", card=None, skip=("cad", "mesh", "fluent", "dimensions"))
        ex = man["parts"]["excel"]
        check(f"(R2-M2) ring_area を 1.01 倍 → 積分の整合 FAIL (相対 {ex['cfd_integral']['rel_to_sum_abs_Q']:.3g})・MANIFEST FAIL・excel が failed_parts に入る",
              not ex["cfd_integral"]["pass"] and man["verdict"] == "FAIL" and "excel" in man["failed_parts"])
        check("(R2-M2) 積分の整合が FAIL なら thermal_use も FAIL", ex["thermal_use"]["verdict"] == "FAIL", ex["thermal_use"].get("reason"))
    # ---- codex result レビュー 3 回目 m1: 単体の CLI が弦の検査を渡し、正常な入力で PASS ----
    import subprocess as _sp
    with tempfile.TemporaryDirectory(prefix="deliv_cli_") as td4:
        cpth = Path(td4) / "card.yaml"; cpth.write_text(yaml.safe_dump(card_ok))
        q = _sp.run([sys.executable, "-m", "forge_design.export.dimension_gate", str(FIX / "step" / "wall_physical.step"),
                     "--card", str(cpth), "--msh", str(FIX / "nozzle.msh")], capture_output=True, text=True,
                    cwd=str(Path(__file__).resolve().parents[1]))
        check("(R3-m1) 単体の寸法ゲート CLI (--msh あり) が正常な入力で PASS・終了 0", q.returncode == 0 and "VERDICT: PASS" in q.stdout,
              q.stdout.strip().splitlines()[-1] if q.stdout.strip() else q.stderr[-200:])
    print(f"FAIL 件数: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
