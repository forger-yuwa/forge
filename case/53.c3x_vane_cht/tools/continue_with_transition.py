#!/usr/bin/env python3
r"""既存の翼 run を**遷移モデル (turbulence.transition: lm2009) つき**で継続する run を作って回す (C3X / Mark II 共通)。

元 run の config・壁温分布・メッシュをそのまま使い、変えるのは (1) `transition: "lm2009"`、(2) 初期場の $k$/$\omega$ (入口値に戻す)、
(3) 任意で入口の $\omega$ (入口粘性比の感度) だけ。plan turbulence-transition-lm2009 §5.1 #6「同じ実測壁温・同じメッシュで遷移 ON/OFF」。

- `--reset-turb` (既定 ON): 流れ場は引き継いで $k$/$\omega$ だけ入口値に戻す。平板 (`case/57`) では戻しても戻さなくても同じ解だったが、
  翼では入口 $\omega$ を変える感度 run で古い乱流場を持ち込まないために戻す。
- `--omega-scale s`: 入口 $\omega$ を s 倍する。入口粘性比 $\mu_t/\mu=\rho k/(\omega\mu)$ は 1/s 倍になる (元 run は約 10)。
- `--off`: 遷移モデルを入れない対照 (同じ出発場・同じ step 数。初期場の違いを切り分けるため)。

usage: continue_with_transition.py SRC_RUN NEW_RUN [--steps 60000] [--out-int 5000] [--omega-scale 1.0] [--no-reset-turb] [--off] [--kato 0|1] [--cfl C] [--dry]
"""
import argparse, os, re, shutil, subprocess, sys
from pathlib import Path
import h5py, numpy as np

ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "solver_density_cuda" / "tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""),
           FORGE_CUDA_BLOCKSIZE="256", FORGE_BIN=str(ROOT / "solver_density_cuda/.build-native/relwithdebinfo/forge"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("src"); ap.add_argument("new")
    ap.add_argument("--steps", type=int, default=60000); ap.add_argument("--out-int", type=int, default=5000)
    ap.add_argument("--omega-scale", type=float, default=1.0); ap.add_argument("--no-reset-turb", action="store_true")
    ap.add_argument("--off", action="store_true"); ap.add_argument("--cfl", type=float, default=None); ap.add_argument("--dry", action="store_true")
    ap.add_argument("--level", type=int, default=None, help="output.level を上書き (2 で遷移モデルの診断場 lm* と残差場も出す)")
    ap.add_argument("--kato", type=int, default=None, help="katoLaunder を上書き (前縁よどみ点の k 過大生成を抑える。LM2009 の著者が併用を勧める)")
    a = ap.parse_args()
    src, new = Path(a.src).resolve(), Path(a.new).resolve()
    if new.exists(): raise SystemExit(f"{new} exists")
    new.mkdir()
    for f in list(src.glob("*.yaml")) + list(src.glob("wall_profile_*.csv")) + list(src.glob("wall_nodes_ordered.csv")) + [src / "mesh.h5"]:
        if f.exists(): shutil.copy(f, new / f.name)
    cfg = (new / "solverConfig.yaml").read_text()
    cfg = re.sub(r"nStepOuter:\s*\d+", f"nStepOuter: {a.steps}", cfg); cfg = re.sub(r"outStepInterval:\s*\d+", f"outStepInterval: {a.out_int}", cfg)
    if a.cfl is not None:
        cfg = re.sub(r"cfl:\s*[0-9.eE+-]+", f"cfl: {a.cfl}", cfg); cfg = re.sub(r"cfl_pseudo:\s*[0-9.eE+-]+", f"cfl_pseudo: {a.cfl}", cfg)
    if not a.off and "transition" not in cfg:     # 継続元が既に遷移つきならそのまま
        cfg, n = re.subn(r'(turbulence:\s*\{[^}]*)\}', r'\1, transition: "lm2009"}', cfg, count=1); assert n == 1
    if a.level is not None:
        cfg, n = re.subn(r"(output:\s*\{[^}]*?level:\s*)\d", rf"\g<1>{a.level}", cfg); assert n == 1
    if a.kato is not None:
        cfg, n = re.subn(r"katoLaunder:\s*\d", f"katoLaunder: {a.kato}", cfg); assert n == 1
    (new / "solverConfig.yaml").write_text(cfg)
    bc = (new / "bcondConfig.yaml").read_text()
    m = re.search(r"\bk:\s*([0-9.eE+-]+),\s*omega:\s*([0-9.eE+-]+)", bc); k_in, om_in = float(m.group(1)), float(m.group(2)) * a.omega_scale
    bc = bc[:m.start()] + f"k: {k_in}, omega: {om_in:.6g}" + bc[m.end():]
    (new / "bcondConfig.yaml").write_text(bc)
    res = sorted(src.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
    r = subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res), str(new / "mesh.h5")], env=ENV, capture_output=True, text=True)
    if r.returncode != 0: raise SystemExit(r.stdout + r.stderr)
    note = f"{res.relative_to(ROOT)}"
    if not a.no_reset_turb:
        with h5py.File(new / "mesh.h5", "r+") as f:
            ro = f["/VALUE/ro"][:].astype(float)
            f["/VALUE/roK"][:] = (ro * k_in).astype(np.float32); f["/VALUE/roOmega"][:] = (ro * om_in).astype(np.float32)
            for n_ in ("roGamma", "roReth"):
                if "/VALUE/" + n_ in f: del f["/VALUE/" + n_]
        note += f"  (k/omega reset to inlet values k={k_in}, omega={om_in:.6g})"
    (new / "CONTINUED_FROM").write_text(note + "\n")
    (new / "GEN_ARGS").write_text(" ".join(sys.argv[1:]) + "\n")
    print("prepared", new.relative_to(ROOT), "|", note)
    if a.dry: return
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(new)], env=ENV, capture_output=True, text=True)
    (new / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    print("rc", r.returncode); print((new / "CONVERGENCE_VERDICT.txt").read_text()[-900:])


if __name__ == "__main__":
    main()
