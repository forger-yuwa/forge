#!/usr/bin/env python3
"""V1 判定: 最終の res_*.h5 で自由流からのずれ (plan §6 V1: 1e-5)。 python3 eval_v1.py RUN_DIR..."""
import glob, re, sys
import h5py, numpy as np, yaml

for run in sys.argv[1:]:
    fs = sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(run + "/res_*.h5") if re.search(r"/res_\d+\.h5$", p))
    st, last = fs[-1]
    fl = yaml.safe_load(open(run + "/bcondConfig.yaml"))["xmin"]["floats"]
    with h5py.File(last, "r") as f:
        V = f["VALUE"]
        ro, P = V["ro"][:].astype(float), V["P"][:].astype(float)
        u = np.stack([V[k][:] for k in ("Ux", "Uy", "Uz")], 1).astype(float)
        uinf = np.array([fl["Ux"], fl["Uy"], fl["Uz"]]); Um = np.linalg.norm(uinf)
        out = [f"dρ {np.max(np.abs(ro / fl['ro'] - 1)):.2e}", f"dP {np.max(np.abs(P / fl['Ps'] - 1)):.2e}",
               f"du {np.max(np.linalg.norm(u - uinf, axis=1)) / Um:.2e}"]
        if "Y0" in V:
            out.append(f"dY {np.max(np.abs(V['Y0'][:] - fl.get('Y0', 0.0))):.2e}")
        if "k" in V and "k" in fl:
            out.append(f"dk {np.max(np.abs(V['k'][:] / fl['k'] - 1)):.2e} (SST は消滅項で一様値は保たれない。判定外)")
        nonfin = sum(int(np.sum(~np.isfinite(V[k][:]))) for k in V.keys())
    worst = max(float(x.split()[1]) for x in out[:3] + ([out[3]] if "Y0" in V else []))
    print(f"{run} step {st}: " + "  ".join(out) + f"  非有限 {nonfin}  → {'PASS' if worst <= 1e-5 and nonfin == 0 else 'FAIL'} (≤ 1e-5)")
