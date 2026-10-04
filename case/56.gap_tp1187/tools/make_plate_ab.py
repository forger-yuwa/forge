#!/usr/bin/env python3
"""2D 平板の壁法線格子 A/B (plan §5.1 #60、codex diagnose 2026-09-26 2 回目) の run を組む。

目的: 3D run の平板部が 2D 平板より 8.9 % 低い (`qwall`、座標を揃えた比較) のが、
**壁法線方向の格子の粗さだけ**で出るかを見る。両腕の違いはメッシュの y 節点列だけにする:

    細側  mesh/fp_yfine.h5   = 既存 fp_t8 と節点集合が一致 (137 節点/列)
    粗側  mesh/fp_ycoarse.h5 = y<0.08 m は 3D (gap3d_c50b_v2) の入口列、その上は細側と同じ (63 節点/列)

段階起動はしない。既存の収束場 (`run_0002_fp_t8_long` の最終 res) から始め、
**本段と同じ設定の 1 段**だけを回す。細側は同一メッシュなので `restart_field.py`、
粗側は `interp_field.py` (cross-mesh)。2D 平板の既存 run と違い、3D run と同じく
`qAccumulatorFP64: 1` を入れる (両腕で揃える。codex Major)。

    python3 tools/make_plate_ab.py --run run_0067_fp_ab_fine   --mesh fp_yfine
    python3 tools/make_plate_ab.py --run run_0068_fp_ab_coarse --mesh fp_ycoarse
    FORGE_CUDA_BLOCKSIZE=128 solver_density_cuda/tools/run_case.sh case/56.gap_tp1187/<run>
"""
import argparse, json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
TOOLS = CASE.parents[1] / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS))
from stage_manifest import StageManifest                         # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", required=True, choices=("fp_yfine", "fp_ycoarse"))
    ap.add_argument("--src", default="run_0002_fp_t8_long/res_150000.h5")
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--out-int", type=int, default=2000)
    a = ap.parse_args()

    rd = CASE / a.run
    if rd.exists() and any(rd.iterdir()):
        raise SystemExit(f"{rd} は既にある。消さない — 別の run 名を使うこと")
    # 入口・壁・物性・IC パッチは既存の 2D 平板と同じ生成器に任せる (--dry = 設定を書くだけ)
    subprocess.run([sys.executable, str(HERE / "make_case.py"), "--run", a.run, "--mesh", a.mesh,
                    "--dry"], check=True)
    st = {s["tag"]: s for s in json.loads((rd / "stages.json").read_text())}["main"]
    st = {**st, "nsteps": a.steps, "out_int": a.out_int}
    import make_case as mc
    cfg = mc.CFG.format(species=", ".join(mc.SPECIES), **st)
    cfg = cfg.replace("detectNaN: 1", "detectNaN: 1, qAccumulatorFP64: 1")
    (rd / "solverConfig.yaml").write_text(cfg, encoding="utf-8")
    (rd / "stages.json").write_text(json.dumps([{**st, "tag": "main"}], indent=2), encoding="utf-8")
    sm = StageManifest(rd)
    sm.add("main", cfg, (rd / "bcondConfig.yaml").read_text(), history="residual_history.csv",
           restart_from=a.src)
    sm.write()

    src = CASE / a.src
    tool = "restart_field.py" if a.mesh == "fp_yfine" else "interp_field.py"
    subprocess.run([sys.executable, str(TOOLS / tool), str(src), str(rd / "mesh.h5")], check=True)
    (rd / "case_ab.json").write_text(json.dumps(dict(
        mesh=a.mesh, init_from=a.src, init_tool=tool, steps=a.steps, out_int=a.out_int,
        cuda_blocksize=128, qAccumulatorFP64=1, plan="§5.1 #60 平板の壁法線格子 A/B"),
        indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{a.run}: {a.mesh}、初期値 {a.src} ({tool})、{a.steps} step")


if __name__ == "__main__":
    main()
