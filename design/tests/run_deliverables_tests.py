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
    print(f"FAIL 件数: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
