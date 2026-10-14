"""plan tooling-nozzle-core-grid §4.104 P-hs: 停止規則の帯 y_b(x) を共通の起点 (run_0916_hi_h_cont の res_40000 = H の本段の通算 10 万 step) で 1 回だけ固定する。
`stop_rule.fix_band` (719dfbf4) をそのまま使う。Euler の参照は M6 の `run_0174_euler_v5d_M_r1` の res_54000・prepare_info を
`~/forge-coregrid/_ext/euler_ref_0174` に複製し、nozzle.h5 は res の MESH から作ったもの (M6 の作業ディレクトリでは格子ファイルが消されていた)。
H の縮流部は Euler の参照 (今の形 C) と形が違うが、停止規則が読むのは x/x_E 0.4・0.7・0.95 (x 15.9〜37.8 r_t、x ≥ 0 の壁は同じ) の断面だけ。
usage: python3 cg_hs_band.py → _band_ab/core_grid/hs_band.json"""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path.home() / "forge-coregrid" / "_ext" / "f5design_719dfbf4" / "design"))
from forge_design.evaluate import stop_rule as SR  # noqa: E402
SRC_RUN, SRC_STEP = "run_0916_hi_h_cont", 40000
EU = Path.home() / "forge-coregrid" / "_ext" / "euler_ref_0174"
out = HERE / "_band_ab" / "core_grid" / "hs_band"
band = SR.fix_band(HERE / SRC_RUN, SRC_STEP, str(EU), out)
band.update(src_run=SRC_RUN, plan="plans/active/tooling-nozzle-core-grid.md §4.104", stop_rule_commit="719dfbf4")
(HERE / "_band_ab" / "core_grid" / "hs_band.json").write_text(json.dumps(band, ensure_ascii=False))
print("band", len(band["x"]), "点", band["seconds"], "s", band["res_sha256"][:16])
