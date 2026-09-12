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
        sf, df = os.path.join(a.src, f), os.path.join(a.dst, f)
        if os.path.exists(sf): shutil.copy(sf, df)
        elif os.path.exists(df): print("using existing", df)          # 元 run にメッシュが無い (別 run から借用) 場合
        else: raise SystemExit("mesh/value file not found: %s (copy it into %s first)" % (sf, a.dst))
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

# 場の差の判定 (plan performance-3d-node-sst-speedup §4.3, 2026-09-12 改訂):
#   正規化: 速度成分 (Ux,Uy,Uz,roUx,roUy,roUz) は速度ベクトルの尺度 max|U| (roU* は max|ρU| ベクトル) で、他は各場の max|a| で割る。
#   絶対基準: A (ro,P,T,roY*,Y*,sonic,vis_lam) ≤1e-5 / B (速度,roe,h0,k,omega,roK,roOmega,g_*,rog_*,roQ*) ≤1e-4 / C (vis_turb) ≤1e-2。
#   wall_dist・res_*・dt_local・限定子など診断量は除外。
CAT_A = {"ro","P","T","sonic","vis_lam"}; CAT_C = {"vis_turb"}
SKIP_PREFIX = ("wall_dist","res_","dt_local","limiter_","ducros","cfl","volume","transport_diag","src_jac","delta_les","l_des","rd_des","fd_","fe_","Pk_diag","Taw_diag","wf_","roK_wf","axisym_divU","condDrdt","condR30","condJ","condTsat","condS_")
import re as _re
GRAD_RE = _re.compile(r"^d[A-Za-z]+d[xyz]$")   # 勾配診断 (dPdx, drody, dUxdz ...) は判定対象外
def category(k):
    if k.startswith(SKIP_PREFIX) or GRAD_RE.match(k): return None
    if k in CAT_A or k.startswith("roY") or (len(k)==2 and k[0]=="Y" and k[1].isdigit()): return "A"
    if k in CAT_C: return "C"
    return "B"
TOL = {"A":1e-5, "B":1e-4, "C":1e-2}
def scales(v):
    import numpy as np
    U = np.sqrt(sum(v[c][...].astype(np.float64)**2 for c in ("Ux","Uy","Uz") if c in v)); rU = np.sqrt(sum(v[c][...].astype(np.float64)**2 for c in ("roUx","roUy","roUz") if c in v))
    return {"Ux":U.max(),"Uy":U.max(),"Uz":U.max(),"roUx":rU.max(),"roUy":rU.max(),"roUz":rU.max()}

def cmp_fields(ref, f, verbose=True):
    """ref/f: h5 VALUE group。判定表を返す (list of (key, cat, normdiff, tol, ok))"""
    sc = scales(ref); rows = []
    for k in sorted(ref.keys()):
        cat = category(k)
        if cat is None: continue
        if k not in f:                       # 比較先に無い場は明示的に FAIL (黙って除外しない, codex result-3 m4)
            rows.append((k, cat, float("inf"), float("inf"), TOL[cat], False, -1)); continue
        x = ref[k][...].astype(np.float64); y = f[k][...].astype(np.float64)
        if x.shape != y.shape:
            rows.append((k, cat, float("inf"), float("inf"), TOL[cat], False, -2)); continue
        s = max(sc.get(k, np.abs(x).max()), 1e-30); d = np.abs(x - y)
        nonfinite = int((~np.isfinite(y)).sum())
        mx = float(d.max()) / s if nonfinite == 0 else float("inf")
        rows.append((k, cat, mx, (np.sqrt((d**2).mean())/s if nonfinite == 0 else float("inf")), TOL[cat], (mx <= TOL[cat]) and nonfinite == 0, nonfinite))
    return rows

def cmp(a):
    """判定: 各場について max 差 ≤ 絶対基準 (TOL) **または** ≤ 2×ノイズ床 (基準バイナリ同士 [ref と --noise ラベル群] の最大差)。
    ノイズ床が絶対基準を超える場 (cell の atomicAdd 床、減衰乱流の速度など) は絶対基準では判定できないので、ノイズ比で判定し
    その旨を表示する。--noise を与えない場合は絶対基準のみ。
    --truth <label>: double ビルド (flowFormat.hpp の typedef を double にした基準コミット) の run ラベル。絶対基準・ノイズ床の両方で
    落ちた場について、**新バイナリの double 解からの距離が基準バイナリの距離以下なら合格** (差が float 丸めの系統差で、基準自身が
    同じだけ double 解からずれている場合。2026-09-12 case/44 run_0200: 等温壁 BL の T は float 版がどれも double 解から 3.1〜3.7e-5)。"""
    cfg = read_cfg(a.dst); n = int(re.search(r"nStepOuter\s*:\s*(\d+)", cfg).group(1))
    ref = h5py.File(os.path.join(a.dst, a.ref, "res_%d.h5" % n), "r")["VALUE"]
    truth_ref = {}
    if a.truth:
        truth = h5py.File(os.path.join(a.dst, a.truth, "res_%d.h5" % n), "r")["VALUE"]
        truth_ref = {k: mx for k, cat, mx, rms, tol, ok, nn in cmp_fields(truth, ref)}
    noise = {}
    for nl in (a.noise or []):
        for k, cat, mx, rms, tol, ok, nn in cmp_fields(ref, h5py.File(os.path.join(a.dst, nl, "res_%d.h5" % n), "r")["VALUE"]):
            noise[k] = max(noise.get(k, 0.0), mx)
    anyfail = False
    for label in a.labels:
        f = h5py.File(os.path.join(a.dst, label, "res_%d.h5" % n), "r")["VALUE"]
        rows = cmp_fields(ref, f); out = []; nbad = 0
        truth_lab = {k: mx for k, cat, mx, rms, tol, ok, nn in cmp_fields(truth, f)} if a.truth else {}
        for k, cat, mx, rms, tol, ok, nn in rows:
            nz = noise.get(k); by_noise = (nz is not None and np.isfinite(mx) and mx <= 2.0 * nz)
            dt_ref, dt_lab = truth_ref.get(k), truth_lab.get(k)
            by_truth = (dt_ref is not None and dt_lab is not None and np.isfinite(dt_lab) and dt_lab <= dt_ref)
            passed = ok or by_noise or by_truth
            if not passed: nbad += 1
            tag = "ok" if ok else ("ok(noise x%.1f)" % (mx / nz) if by_noise else ("ok(truth: new %.2e <= base %.2e)" % (dt_lab, dt_ref) if by_truth else ("MISSING" if nn == -1 else "SHAPE" if nn == -2 else ("NONFINITE(%d)" % nn) if nn > 0 else "EXCEED")))
            if a.truth and dt_ref is not None and not (ok or by_noise): out.append("  %-10s      truth-dist base=%.2e new=%.2e" % (k, dt_ref, dt_lab if dt_lab is not None else float("nan")))
            out.append("  %-10s [%s] max=%.2e rms=%.2e tol=%.0e%s %s" % (k, cat, mx, rms, tol, (" noise=%.2e" % nz) if nz is not None else "", tag))
        print("--- %s vs %s : %s (%d/%d fields; noise labels: %s%s)" % (a.ref, label, "PASS" if nbad == 0 else "FAIL", len(rows) - nbad, len(rows), ",".join(a.noise or []), ("; truth: " + a.truth) if a.truth else ""))
        for o in out: print(o)
        anyfail = anyfail or (nbad > 0)
    if anyfail: sys.exit(1)

ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
p = sp.add_parser("prep"); p.add_argument("src"); p.add_argument("dst"); p.add_argument("ic"); p.add_argument("--mesh"); p.add_argument("--value"); p.add_argument("--nsteps", type=int, default=300); p.set_defaults(fn=prep)
p = sp.add_parser("run"); p.add_argument("dst"); p.add_argument("--bin", action="append", required=True); p.add_argument("--cfgsub"); p.set_defaults(fn=run)
p = sp.add_parser("cmp"); p.add_argument("dst"); p.add_argument("ref"); p.add_argument("labels", nargs="+"); p.add_argument("--noise", action="append", help="基準バイナリの別 run ラベル (ノイズ床用, 複数可)"); p.add_argument("--truth", help="double ビルドの run ラベル (絶対/ノイズで落ちた場は double 解からの距離が基準以下なら合格)"); p.set_defaults(fn=cmp)
a = ap.parse_args(); a.fn(a)
