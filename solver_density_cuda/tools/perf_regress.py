#!/usr/bin/env python3
"""既存 run の場から短時間継続して 2 バイナリ (base/new) の場を比較する回帰ドライバ
(plan performance-3d-node-sst-speedup §6)。

  perf_regress.py prep SRC_RUN DST_RUN IC_RES [--mesh NAME] [--value NAME] [--nsteps N]
  perf_regress.py run  DST_RUN --bin LABEL=PATH [--bin LABEL=PATH ...] [--cfgsub "old=new"]
  perf_regress.py cmp  DST_RUN REF_LABEL LABEL [LABEL ...]

prep: SRC_RUN の config/bcond/species/probe とメッシュ h5 (valueFileName) を DST_RUN に複製し、IC_RES の VALUE を index コピー
      (保存量 ro* が無ければ作る)。nStepOuter=N, outStepInterval=N に書き換える。
run : DST_RUN/<LABEL>/ に forge を回す (res_N.h5 が残る)。
cmp : 場のスケール max|a| で正規化した最大差を全 VALUE について表示。"""
import argparse, os, re, shutil, subprocess, sys
import h5py, numpy as np

def read_cfg(d):
    return open(os.path.join(d, "solverConfig.yaml")).read()

def prep(a):
    os.makedirs(a.dst, exist_ok=True)
    cfg = read_cfg(a.src)
    mesh = a.mesh or re.search(r'meshFileName:\s*"?([^"\s]+)"?', cfg).group(1)
    value = a.value or re.search(r'valueFileName:\s*"?([^"\s]+)"?', cfg).group(1)
    for f in ["bcondConfig.yaml", "species_db.yaml", "probe.yaml"]:
        if os.path.exists(os.path.join(a.src, f)): shutil.copy(os.path.join(a.src, f), a.dst)
    for f in {mesh, value}:
        shutil.copy(os.path.join(a.src, f), os.path.join(a.dst, f))
    with h5py.File(a.ic, "r") as s, h5py.File(os.path.join(a.dst, value), "r+") as d:
        n = len(d["VALUE/ro"]); nk = 0; created = []
        for k in s["VALUE"]:
            if k == "wall_dist" or len(s["VALUE"][k]) != n: continue
            if k not in d["VALUE"]:
                if not k.startswith("ro"): continue
                d["VALUE"].create_dataset(k, data=s["VALUE"][k][...]); created.append(k); nk += 1; continue
            d["VALUE"][k][:] = s["VALUE"][k][...]; nk += 1
    cfg = re.sub(r"(nStepOuter\s*:\s*)\d+", r"\g<1>%d" % a.nsteps, cfg, count=1)
    cfg = re.sub(r"(outStepInterval\s*:\s*)\d+", r"\g<1>%d" % a.nsteps, cfg, count=1)
    open(os.path.join(a.dst, "solverConfig.yaml"), "w").write(cfg)
    open(os.path.join(a.dst, "IC_FROM.txt"), "w").write("%s (index copy, %d fields, created %s)\n" % (os.path.abspath(a.ic), nk, created))
    print("prepared", a.dst, "n", n, "fields", nk, "created", created)

def run(a):
    env = dict(os.environ); env["LD_LIBRARY_PATH"] = "/usr/lib/x86_64-linux-gnu/hdf5/serial:" + env.get("LD_LIBRARY_PATH", "")
    for spec in a.bin:
        label, path = spec.split("=", 1)
        r = os.path.join(a.dst, label); os.makedirs(r, exist_ok=True)
        for f in os.listdir(a.dst):
            p = os.path.join(a.dst, f)
            if not os.path.isfile(p): continue
            if f.endswith(".h5"):
                if os.path.exists(os.path.join(r, f)): os.remove(os.path.join(r, f))
                try: os.link(p, os.path.join(r, f))
                except OSError: shutil.copy(p, os.path.join(r, f))
            elif f.endswith((".yaml", ".txt")): shutil.copy(p, r)
        if a.cfgsub:
            old, new = a.cfgsub.split("=", 1); old = old.encode().decode("unicode_escape"); new = new.encode().decode("unicode_escape")
            t = read_cfg(r); assert old in t, "cfgsub old text not found"; open(os.path.join(r, "solverConfig.yaml"), "w").write(t.replace(old, new))
        with open(os.path.join(r, "forge_run.log"), "w") as log:
            rc = subprocess.call([path], cwd=r, stdout=log, stderr=subprocess.STDOUT, env=env)
        t = open(os.path.join(r, "forge_run.log")).read()
        m = re.search(r"^Time = .*$", t, re.M)
        print("%s/%s exit=%d %s" % (os.path.basename(a.dst), label, rc, m.group(0) if m else "(no Time line)"))

def cmp(a):
    cfg = read_cfg(a.dst); n = int(re.search(r"nStepOuter\s*:\s*(\d+)", cfg).group(1))
    ref = h5py.File(os.path.join(a.dst, a.ref, "res_%d.h5" % n), "r")["VALUE"]
    for label in a.labels:
        f = h5py.File(os.path.join(a.dst, label, "res_%d.h5" % n), "r")["VALUE"]
        print("--- %s vs %s" % (a.ref, label))
        worst = 0.0
        for k in sorted(ref.keys()):
            if k not in f or k == "wall_dist": continue
            x = ref[k][...].astype(np.float64); y = f[k][...].astype(np.float64)
            s = max(np.abs(x).max(), 1e-30); d = np.abs(x - y)
            print("  %-10s max|d|/max|a|=%.2e rms=%.2e nan=%d" % (k, d.max() / s, np.sqrt((d ** 2).mean()) / s, int(np.isnan(y).sum())))

ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
p = sp.add_parser("prep"); p.add_argument("src"); p.add_argument("dst"); p.add_argument("ic"); p.add_argument("--mesh"); p.add_argument("--value"); p.add_argument("--nsteps", type=int, default=300); p.set_defaults(fn=prep)
p = sp.add_parser("run"); p.add_argument("dst"); p.add_argument("--bin", action="append", required=True); p.add_argument("--cfgsub"); p.set_defaults(fn=run)
p = sp.add_parser("cmp"); p.add_argument("dst"); p.add_argument("ref"); p.add_argument("labels", nargs="+"); p.set_defaults(fn=cmp)
a = ap.parse_args(); a.fn(a)
