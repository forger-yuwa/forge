"""成果物の一括書き出し `forge_design.export.deliverables` の試験 (plan tooling-nozzle-wt-deliverables §6 V6・V7 のゲート部分)。

fixture は FORGE_DELIV_FIXTURE (既定: このセッションの scratchpad の m6_prod; 作り方は tests/make_deliverables_fixture.py)。
M6 の生産の壁の準備ディレクトリ (wall_repr.json・prepare_info.json・bcondConfig.yaml・nozzle.msh) を仮の run にして、
(a) 一括書き出しが通り MANIFEST の sha256 が実ファイルと一致、(b) 寸法ゲート: 生産の設計値で作ったカードは PASS・出口長を 2.2 mm
ずらしたカードは FAIL で MANIFEST の verdict も FAIL、(c) 入力が欠けた run (wall_repr.json なし) は cad を理由つきで飛ばす、を確かめる。
FreeCAD・gmsh が要る (無ければ飛ばす)。作業は一時ディレクトリで行い、fixture には書かない。
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


def main() -> int:
    if not (FIX / "wall_repr.json").exists() or not (FIX / "nozzle.msh").exists():
        print(f"skip: fixture が無い ({FIX})。tests/make_deliverables_fixture.py で作る")
        return 0
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
    print(f"FAIL 件数: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
