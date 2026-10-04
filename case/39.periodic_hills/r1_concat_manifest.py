#!/usr/bin/env python3
"""R1 延長 run に「元 run の S1_spin → S2_dev + 延長 S3_ext」の stage_manifest を書く
(plan boundary-node-periodic-gradient-fix §5.1 #6d、codex result M5)。

延長 run (run_0038 / run_0039) は元 run (run_0036 / run_0037) の res_800000 を `restart_field.py` で写して
同一設定で継続したもの (差は valueFileName とコメントだけ)。`--from-floor` の判定では延長開始直後の
スパイクだけで NOT CONVERGED になったので、**同一実効設定の 1 区間**として連結し `check_convergence.py --segment`
で判定する (スパイクを含む全履歴が判定区間に入る)。区間は名前でなく hard キー (方程式・BC・離散化) で決まる
(`stage_manifest.py`)。元 run の段の履歴は相対パス (`../<元 run>/residual_history_<tag>.csv`) で参照し、コピーしない。

    python3 r1_concat_manifest.py EXT_RUN ORIG_RUN
    python3 ../../solver_density_cuda/tools/stage_manifest.py EXT_RUN --segments
    python3 ../../solver_density_cuda/tools/check_convergence.py EXT_RUN --segment
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "solver_density_cuda", "tools"))
from stage_manifest import StageManifest  # noqa: E402


def main():
    ext, orig = (os.path.abspath(p) for p in sys.argv[1:3])
    rel = os.path.relpath(orig, ext)
    if os.path.exists(os.path.join(ext, "stage_manifest.json")):
        raise SystemExit(f"既に stage_manifest.json がある: {ext}")
    sm = StageManifest(ext)
    bc_orig = open(os.path.join(orig, "bcondConfig.yaml")).read()
    for tag in ("S1_spin", "S2_dev"):
        cfg = open(os.path.join(orig, f"solverConfig_{tag}.yaml")).read()
        sm.add(tag, cfg, bc_orig, history=os.path.join(rel, f"residual_history_{tag}.csv"))
    sm.add("S3_ext", open(os.path.join(ext, "solverConfig.yaml")).read(),
           open(os.path.join(ext, "bcondConfig.yaml")).read(), history="residual_history.csv", restart_from="S2_dev")
    print("wrote", sm.write())


if __name__ == "__main__":
    main()
