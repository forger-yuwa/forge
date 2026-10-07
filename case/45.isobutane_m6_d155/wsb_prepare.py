"""plan tooling-nozzle-wall-single-bspline §5.1 #2b・§6 W0/W3: 固定入力で `prepare_ns` まで回し (forge は起動しない)、
成果物のハッシュを記録する。

固定入力 = 単調壁の生産問題 `problem_d155_ns_finemesh_recal_final_mono.yaml` (δ_r は YAML の `deltastar_initializer` の積分法経路 =
run_0147 と同じ)。`initial_line_run` だけを凍結源の絶対パス (CASE_RUNS 以下) に書き換えた写しを <out>/../inputs/ に置いて渡す。
prepare_ns の引数は prep_c2pin.py と同じ (nsteps 12000・cfl_main 5.0・implicit_relax 0.7、IC なし)。

usage:
  [DESIGN_DIR=<design/ の写し>] [CASE_RUNS=<run_0062・run_0147 のある case dir>] \
  FORGE_CONVERTER=<変換器> FORGE_ALLOW_UNVERIFIED_SPECIES=1 \
  python wsb_prepare.py OUT_DIR [--repr legacy|single_bspline|none] [--label TEXT]
    --repr none (既定): YAML にキーを書かない (変更前のコードの基準・W0)。legacy / single_bspline: geometry.physical_wall_repr を書く。
  → OUT_DIR (prepare_ns の成果物) と OUT_DIR/../<OUT_DIR 名>.hashes.json (ファイルの sha256・HDF5 のデータセット単位の sha256・環境)

HDF5 はファイルのバイト列が同じ入力でも再現しない (変換器の既知の性質) ので、比較はデータセット単位で行う (wsb_compare.py)。
"""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
DESIGN = Path(os.environ.get("DESIGN_DIR", str(C.parents[1] / "design"))).resolve()
sys.path.insert(0, str(DESIGN))
RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
PN = "problem_d155_ns_finemesh_recal_final_mono.yaml"
REF_RUN = RUNS / "run_0147_ns_mono_final"


def sha_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def h5_digest(path: Path) -> dict:
    """HDF5 のデータセットごとの sha256 (dtype・shape・データのバイト列) と属性 (グループ・データセット)。"""
    import h5py
    out = {"datasets": {}, "attrs": {}}

    def _attrs(name, obj):
        if len(obj.attrs):
            out["attrs"][name or "/"] = {k: (v.tolist() if hasattr(v, "tolist") else (v.decode() if isinstance(v, bytes) else v))
                                         for k, v in obj.attrs.items()}

    with h5py.File(path, "r") as f:
        _attrs("/", f)

        def visit(name, obj):
            _attrs(name, obj)
            if isinstance(obj, h5py.Dataset):
                a = np.asarray(obj[()])
                h = hashlib.sha256()
                h.update(str(a.dtype).encode()); h.update(str(a.shape).encode())
                h.update(np.ascontiguousarray(a).tobytes() if a.dtype.kind != "O" else repr(a.tolist()).encode())
                out["datasets"][name] = {"dtype": str(a.dtype), "shape": list(a.shape), "sha256": h.hexdigest()}
        f.visititems(visit)
    return out


def problem_copy(out_dir: Path, repr_key: str) -> Path:
    """固定入力の問題 YAML の写し (initial_line_run を絶対パスに; repr_key != none なら geometry.physical_wall_repr を足す)。"""
    txt = (C / PN).read_text()
    il = "initial_line_run: run_0062_euler_wallfit_fit_r1_ext6k"
    if txt.count(il) != 1:
        raise RuntimeError(f"{PN} に {il!r} が 1 回だけ無い")
    txt = txt.replace(il, f"initial_line_run: {RUNS / 'run_0062_euler_wallfit_fit_r1_ext6k'}")
    if repr_key != "none":
        anchor = "  pw_ramp: [-11.0, -6.0]"
        if txt.count(anchor) != 1:
            raise RuntimeError("pw_ramp の行が 1 回だけ無い")
        txt = txt.replace(anchor, anchor + f"\n  physical_wall_repr: {repr_key}", 1)
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
    ap.add_argument("--repr", default="none", choices=("none", "legacy", "single_bspline"))
    ap.add_argument("--label", default="")
    a = ap.parse_args(argv)
    out_dir = Path(a.out_dir).resolve()
    if out_dir.exists():
        raise SystemExit(f"{out_dir} が既にある (上書きしない)")
    conv = os.environ.get("FORGE_CONVERTER")
    if not conv or not Path(conv).is_file():
        raise SystemExit("FORGE_CONVERTER (変換器) を指定する")
    prob = problem_copy(out_dir, a.repr)
    import scipy
    from forge_design.evaluate.runner_axismach import prepare_ns
    info = prepare_ns(prob, out_dir, nsteps=12000, ic_from=None, cfl_main=5.0, implicit_relax=0.7)
    files = {}
    for f in sorted(out_dir.iterdir()):
        if f.is_file():
            files[f.name] = sha_file(f)
    rec = {
        "label": a.label, "out_dir": str(out_dir), "repr_key": a.repr, "problem_copy": str(prob), "problem_copy_sha256": sha_file(prob),
        "problem_source": str(C / PN), "problem_source_sha256": sha_file(C / PN),
        "design_dir": str(DESIGN), "python": sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
        "converter": conv, "converter_sha256": sha_file(Path(conv)),
        "FORGE_ALLOW_UNVERIFIED_SPECIES": os.environ.get("FORGE_ALLOW_UNVERIFIED_SPECIES"),
        "initial_line_run": str(RUNS / "run_0062_euler_wallfit_fit_r1_ext6k"),
        "files_sha256": files,
        "h5": {"nozzle.h5": h5_digest(out_dir / "nozzle.h5")},
        "prepare_info_keys": sorted(info),
    }
    try:
        rec["design_git"] = subprocess.run(["git", "-C", str(DESIGN), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
    except OSError:
        rec["design_git"] = None
    # run_0147 の保存した δ_r 表との一致 (prepare_ns は initializer 経路で δ_r を作り直す)
    ref = REF_RUN / "delta_r_initial.csv"
    new = out_dir / "delta_r_initial.csv"
    if ref.is_file() and new.is_file():
        ta = np.loadtxt(ref, delimiter=",", skiprows=1); tb = np.loadtxt(new, delimiter=",", skiprows=1)
        rec["delta_r_initial_vs_run_0147"] = {
            "ref": str(ref), "ref_sha256": sha_file(ref), "file_identical": sha_file(ref) == sha_file(new),
            "shape": [list(ta.shape), list(tb.shape)],
            "max_abs_diff_all_cols": (float(np.abs(ta - tb).max()) if ta.shape == tb.shape else None),
            "max_abs_diff_delta_r_rt": (float(np.abs(ta[:, 1] - tb[:, 1]).max()) if ta.shape == tb.shape else None),
            "x_identical": bool(ta.shape == tb.shape and np.array_equal(ta[:, 0], tb[:, 0]))}
    else:
        rec["delta_r_initial_vs_run_0147"] = {"ref_exists": ref.is_file(), "new_exists": new.is_file()}
    wp = REF_RUN / "wall_physical.csv"
    if wp.is_file() and (out_dir / "wall_physical.csv").is_file():
        ta = np.loadtxt(wp, delimiter=",", skiprows=1); tb = np.loadtxt(out_dir / "wall_physical.csv", delimiter=",", skiprows=1)
        rec["wall_physical_vs_run_0147"] = {"file_identical": sha_file(wp) == sha_file(out_dir / "wall_physical.csv"),
                                            "max_abs_diff_m": (float(np.abs(ta - tb).max()) if ta.shape == tb.shape else None)}
    # W4 統合検査に使う既存結果 (run_0147 の nozzle.h5・res_*.h5) はローカルに無い (AWS)。取得元を記録するだけで取りに行かない
    rec["run_0147_results"] = {"local_dir": str(REF_RUN), "nozzle.h5_local": (REF_RUN / "nozzle.h5").is_file(),
                               "res_h5_local": sorted(p.name for p in REF_RUN.glob("res_*.h5")),
                               "source": ("AWS 共有 GPU インスタンス上の case/45.isobutane_m6_d155/run_0147_ns_mono_final/ "
                                          "(run_mono_ns_chain.sh で作成、延長は run_0149_ns_mono_final_ext)。未取得 — "
                                          "リモートの正確なパス・ファイル名・sha256 は取りに行っていない")}
    dst = out_dir.parent / f"{out_dir.name}.hashes.json"
    dst.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({k: rec[k] for k in ("out_dir", "repr_key", "files_sha256", "delta_r_initial_vs_run_0147")}, indent=1, default=str))


if __name__ == "__main__":
    main()
