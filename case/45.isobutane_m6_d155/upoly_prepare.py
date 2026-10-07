"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U0・U3 (と全域 1 本の plan の W3 のやり直し): 固定入力で `prepare_ns` まで回し
(forge は起動しない)、成果物のハッシュを記録する。wsb_prepare.py (全域 1 本の plan) の流用。

固定入力 = 単調壁の生産問題 `problem_d155_ns_finemesh_recal_final_mono.yaml` (δ_r は YAML の `deltastar_initializer` の積分法経路 =
run_0147 と同じ)。`initial_line_run` を凍結源の絶対パス (CASE_RUNS 以下) に書き換えた写しを <out>/../inputs/ に置いて渡す。
prepare_ns の引数は prep_c2pin.py と同じ (nsteps 12000・cfl_main 5.0・implicit_relax 0.7、IC なし)。

variant (問題の写しの書き換え):
  ramp_explicit       : そのまま (移行で pw_upstream: ramp + pw_ramp [-11, -6] を明示済み) — U0 の変更後
  poly_default        : pw_ramp・pw_upstream の行を消す (キー無し = 既定 poly)
  poly_legacy         : pw_ramp を消し pw_upstream: poly と physical_wall_repr: legacy を書く — W3 の腕 A
  poly_single_bspline : pw_ramp を消し pw_upstream: poly と physical_wall_repr: single_bspline を書く — W3 の腕 B
  任意で --r-throat R (spec.r_throat を書き換え) と --sizing METHOD:TARGET_M (spec.sizing を書く) — U2 の往復

usage:
  [CASE_RUNS=<run_0062 のある case dir>] FORGE_CONVERTER=<変換器> FORGE_ALLOW_UNVERIFIED_SPECIES=1 \
  python upoly_prepare.py OUT_DIR --variant VARIANT [--label TEXT] [--r-throat R] [--sizing throat:0.0767] [--delta-r-csv CSV]
  → OUT_DIR (prepare_ns の成果物) と OUT_DIR/../<OUT_DIR 名>.hashes.json
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import wsb_prepare as WP  # noqa: E402  (DESIGN_DIR を sys.path に入れる・ハッシュの関数)

RUNS = WP.RUNS
PN = WP.PN
VARIANTS = ("ramp_explicit", "poly_default", "poly_legacy", "poly_single_bspline")
RAMP_LINE = re.compile(r"^  pw_ramp: \[-11\.0, -6\.0\].*\n", re.M)
UP_LINE = re.compile(r"^  pw_upstream: ramp .*\n", re.M)


def problem_copy(out_dir: Path, variant: str, r_throat=None, sizing=None) -> Path:
    txt = (C / PN).read_text()
    il = "initial_line_run: run_0062_euler_wallfit_fit_r1_ext6k"
    if txt.count(il) != 1:
        raise RuntimeError(f"{PN} に {il!r} が 1 回だけ無い")
    txt = txt.replace(il, f"initial_line_run: {RUNS / 'run_0062_euler_wallfit_fit_r1_ext6k'}")
    if len(RAMP_LINE.findall(txt)) != 1 or len(UP_LINE.findall(txt)) != 1:
        raise RuntimeError("pw_ramp・pw_upstream: ramp の行が 1 回ずつ無い (移行前の YAML?)")
    if variant != "ramp_explicit":
        txt = UP_LINE.sub("", RAMP_LINE.sub("", txt))
        add = {"poly_default": "", "poly_legacy": "  pw_upstream: poly\n  physical_wall_repr: legacy\n",
               "poly_single_bspline": "  pw_upstream: poly\n  physical_wall_repr: single_bspline\n"}[variant]
        anchor = "  wall_fit_mono_r2: [0.0, 1.5]"
        i = txt.index(anchor)
        j = txt.index("\n", i) + 1
        txt = txt[:j] + add + txt[j:]
    if r_throat is not None:
        m = re.findall(r"^  r_throat: [0-9.eE+-]+", txt, re.M)
        if len(m) != 1:
            raise RuntimeError("spec.r_throat の行が 1 回だけ無い")
        txt = re.sub(r"^  r_throat: [0-9.eE+-]+", f"  r_throat: {float(r_throat)!r}", txt, count=1, flags=re.M)
    if sizing is not None:
        meth, tgt = sizing.split(":")
        txt = re.sub(r"^(  r_throat: .*\n)", lambda mm: mm.group(1) + f"  sizing: {{method: {meth}, target_m: {float(tgt)!r}}}\n",
                     txt, count=1, flags=re.M)
    inp = out_dir.parent / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    dst = inp / f"{out_dir.name}.problem.yaml"
    if dst.exists():
        raise SystemExit(f"{dst} が既にある (上書きしない)")
    dst.write_text(txt)
    return dst


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--variant", required=True, choices=VARIANTS)
    ap.add_argument("--label", default="")
    ap.add_argument("--r-throat", type=float, default=None)
    ap.add_argument("--sizing", default=None)
    ap.add_argument("--delta-r-csv", default=None, help="prepare_ns(delta_r_csv=...) (NS 後の補正の表から壁を作る)")
    a = ap.parse_args(argv)
    out_dir = Path(a.out_dir).resolve()
    if out_dir.exists():
        raise SystemExit(f"{out_dir} が既にある (上書きしない)")
    conv = os.environ.get("FORGE_CONVERTER")
    if not conv or not Path(conv).is_file():
        raise SystemExit("FORGE_CONVERTER (変換器) を指定する")
    prob = problem_copy(out_dir, a.variant, a.r_throat, a.sizing)
    import scipy
    from forge_design.evaluate.runner_axismach import prepare_ns
    kw = {} if a.delta_r_csv is None else {"delta_r_csv": a.delta_r_csv, "offset": "radial"}
    info = prepare_ns(prob, out_dir, nsteps=12000, ic_from=None, cfl_main=5.0, implicit_relax=0.7, **kw)
    files = {f.name: WP.sha_file(f) for f in sorted(out_dir.iterdir()) if f.is_file()}
    rec = {"label": a.label, "out_dir": str(out_dir), "variant": a.variant, "r_throat": a.r_throat, "sizing": a.sizing,
           "delta_r_csv": a.delta_r_csv, "problem_copy": str(prob), "problem_copy_sha256": WP.sha_file(prob),
           "problem_source": str(C / PN), "problem_source_sha256": WP.sha_file(C / PN),
           "design_dir": str(WP.DESIGN), "python": sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
           "converter": conv, "converter_sha256": WP.sha_file(Path(conv)),
           "FORGE_ALLOW_UNVERIFIED_SPECIES": os.environ.get("FORGE_ALLOW_UNVERIFIED_SPECIES"),
           "initial_line_run": str(RUNS / "run_0062_euler_wallfit_fit_r1_ext6k"), "files_sha256": files,
           "h5": {"nozzle.h5": WP.h5_digest(out_dir / "nozzle.h5")}, "prepare_info_keys": sorted(info),
           "pw_upstream": info.get("pw_upstream"), "sizing_record": info.get("sizing")}
    try:
        rec["design_git"] = subprocess.run(["git", "-C", str(WP.DESIGN), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
        rec["design_git_dirty"] = bool(subprocess.run(["git", "-C", str(WP.DESIGN), "status", "--porcelain", "--", "."],
                                                      capture_output=True, text=True).stdout.strip())
    except OSError:
        rec["design_git"] = None
    dst = out_dir.parent / f"{out_dir.name}.hashes.json"
    dst.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({k: rec[k] for k in ("out_dir", "variant", "pw_upstream", "sizing_record")}, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
