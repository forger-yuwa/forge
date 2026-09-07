#!/usr/bin/env python3
"""相似メッシュのスケール不変性テスト (codex レビュー 2026-09-08 の提案)。

同じ .msh を座標だけ α 倍した相似メッシュで、粘性 μ・熱伝導 λ を α 倍 (Re, Pr 不変)、
その他の物理条件 (ρ, P, T, u, k) を同じにして N step 回すと、**正しくスケール不変な離散化**なら
場 (ρ, P, T, u, k) は一致し ω は 1/α 倍になる。単位付きの閾値・床 (例: 旧 scalarTransport の
|d·S| < 1e-12 [m³]) が混ざっていると α に依存して破綻する。

  usage: test_scale_invariance.py <run_dir_template> [--alpha 1e-3] [--steps 50] [--tol 1e-3]

run_dir_template: solverConfig.yaml / bcondConfig.yaml / (species_db.yaml, probe.yaml) と、
config の meshFileName に対応する **.msh** を含むディレクトリ (変換はこのツールが行う)。
viscMethod は 0 (定数粘性) であること (Sutherland は温度依存でスケールしないため)。
出力: <run_dir_template>_scale_a1 / _scale_aALPHA を作り、比較結果と VERDICT を表示。
"""
import argparse, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np, h5py, yaml

HERE = Path(__file__).resolve().parent
BUILD = HERE.parent / "build"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial")


def scale_msh(src: Path, dst: Path, alpha: float):
    """gmsh msh4.1 の $Nodes 座標だけを α 倍して書き出す。"""
    lines = src.read_text().splitlines()
    out = []; i = 0; n = len(lines)
    while i < n:
        out.append(lines[i])
        if lines[i].strip() == "$Nodes":
            hdr = lines[i + 1].split(); out.append(lines[i + 1]); i += 2
            nblocks = int(hdr[0])
            for _ in range(nblocks):
                bh = lines[i].split(); out.append(lines[i]); i += 1
                nn = int(bh[3])
                for _ in range(nn): out.append(lines[i]); i += 1            # node tags
                for _ in range(nn):
                    x, y, z = (float(v) for v in lines[i].split()[:3]); out.append(f"{x*alpha:.10g} {y*alpha:.10g} {z*alpha:.10g}"); i += 1
            continue
        i += 1
    dst.write_text("\n".join(out) + "\n")


def prepare(tpl: Path, alpha: float, steps: int) -> Path:
    tag = "a1" if alpha == 1.0 else f"a{alpha:g}"
    d = tpl.parent / (tpl.name + f"_scale_{tag}"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
    cfg = yaml.safe_load((tpl / "solverConfig.yaml").read_text())
    h5 = cfg["mesh"]["meshFileName"]; msh = [f for f in tpl.glob("*.msh")][0]
    for f in ["bcondConfig.yaml", "species_db.yaml", "probe.yaml"]:
        if (tpl / f).exists(): shutil.copy(tpl / f, d / f)
    txt = (tpl / "solverConfig.yaml").read_text()
    if cfg["physProp"].get("viscMethod", 0) != 0:
        raise SystemExit("viscMethod must be 0 (constant viscosity) for the similarity test")
    for key in ("visc", "thermCond"):
        v = float(cfg["physProp"][key]); txt = re.sub(rf"(\b{key}:\s*)[0-9.eE+-]+", lambda m: f"{m.group(1)}{v*alpha:.8g}", txt, count=1)
    # 時間の次元を持つ設定 (dt 種, dt_min, dt_max) も α 倍 (局所 dt ∝ L/c なので、これを忘れると dt_min に張り付いて相似が壊れる)
    for key in ("dt", "dt_min", "dt_max"):
        if key in cfg["time"]["deltaT"]:
            v = float(cfg["time"]["deltaT"][key]); txt = re.sub(rf"(\b{key}:\s*)[0-9.eE+-]+", lambda m: f"{m.group(1)}{v*alpha:.8g}", txt, count=1)
    txt = re.sub(r"nStepOuter: \d+", f"nStepOuter: {steps}", txt); txt = re.sub(r"outStepInterval: \d+", f"outStepInterval: {steps}", txt)
    # ω の次元は 1/s: 境界条件 (bcond floats の omega) と config の omegaInf を 1/α 倍する
    txt = re.sub(r"(\bomegaInf:\s*)([0-9.eE+-]+)", lambda m: f"{m.group(1)}{float(m.group(2))/alpha:.8g}", txt)
    (d / "solverConfig.yaml").write_text(txt)
    if (d / "bcondConfig.yaml").exists():
        b = (d / "bcondConfig.yaml").read_text()
        b = re.sub(r"(\bomega:\s*)([0-9.eE+-]+)", lambda m: f"{m.group(1)}{float(m.group(2))/alpha:.8g}", b)
        (d / "bcondConfig.yaml").write_text(b)
    scale_msh(msh, d / msh.name, alpha)
    r = subprocess.run([str(BUILD / "convertGmshToForge"), msh.name, h5], cwd=d, env=ENV, capture_output=True, text=True)
    if not (d / h5).exists(): raise SystemExit("convert failed:\n" + r.stdout[-2000:] + r.stderr[-2000:])
    # IC: テンプレートの h5 に VALUE があればそれを index コピー (ω は 1/α)
    if (tpl / h5).exists():
        with h5py.File(tpl / h5, "r") as s, h5py.File(d / h5, "r+") as t:
            n = len(t["VALUE/ro"])
            for k in s["VALUE"]:
                if k != "wall_dist" and k in t["VALUE"] and len(s["VALUE"][k]) == n:
                    v = s["VALUE"][k][:]
                    if k in ("roOmega", "omega"): v = v / alpha
                    t["VALUE"][k][:] = v
    return d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("template"); ap.add_argument("--alpha", type=float, default=1e-3)
    ap.add_argument("--steps", type=int, default=50); ap.add_argument("--tol", type=float, default=1e-3); a = ap.parse_args()
    tpl = Path(a.template).resolve()
    runs = {alpha: prepare(tpl, alpha, a.steps) for alpha in (1.0, a.alpha)}
    for alpha, d in runs.items():
        subprocess.run([str(HERE / "run_case.sh"), str(d)], env=ENV, capture_output=True, text=True)
    res = {alpha: sorted(d.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1] for alpha, d in runs.items()}
    ok = True
    with h5py.File(res[1.0], "r") as f1, h5py.File(res[a.alpha], "r") as f2:
        for k in ["ro", "P", "T", "Ux", "Uy", "Uz", "k", "omega", "vis_turb"]:
            if f"VALUE/{k}" not in f1: continue
            x = f1[f"VALUE/{k}"][:].astype(float); y = f2[f"VALUE/{k}"][:].astype(float)
            if k == "omega": y = y * a.alpha          # ω ~ 1/L
            if k == "vis_turb": y = y / a.alpha       # μt ~ ρ k/ω ~ L
            den = max(np.max(np.abs(x)), 1e-30); rel = np.max(np.abs(x - y)) / den
            flag = "OK " if rel <= a.tol else "NG "
            if rel > a.tol: ok = False
            j = int(np.argmax(np.abs(x - y))); c1 = f1["MESH/COORD"][:].reshape(-1, 3)[j]
            print(f"  {flag}{k:9s} max|Δ|/max = {rel:.3e}  (at x={c1[0]:.4g} y={c1[1]:.4g}: {x[j]:.5g} vs {y[j]:.5g})")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'} (alpha={a.alpha:g}, steps={a.steps}, tol={a.tol:g})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
