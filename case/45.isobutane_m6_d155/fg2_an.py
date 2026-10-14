"""plan architecture-float-state-double-geometry §6.3 (段 ②) の合格条件 1〜5 の集計。fg2.sh の後に AWS 上で回す。
判定は事前登録 (§6.3) のとおり。結果を標準出力と fg2_verdict.json に書く。"""
import csv, json, h5py, numpy as np
from pathlib import Path
H = Path("/home/ubuntu/forge-wallfit/case")
C45, C56, C39, C48 = H / "45.isobutane_m6_d155", H / "56.gap_tp1187", H / "39.periodic_hills", H / "48.flat_plate_cooled_m4"
DUMPS = {"case/45": (C45 / "run_0396_fg2_dump32", C45 / "run_0397_fg2_dump64"),
         "case/56": (C56 / "run_0082_fg2_dump32", C56 / "run_0083_fg2_dump64"),
         "case/39": (C39 / "run_0956_fg2_dump32", C39 / "run_0957_fg2_dump64"),
         "case/48": (C48 / "run_0046_fg2_dump32", C48 / "run_0047_fg2_dump64")}
V = {}
def verdict(key, ok, detail):
    V[key] = {"pass": bool(ok), "detail": detail}
    print(f"  [{'PASS' if ok else 'FAIL'}] {key}: {detail}")
def ds(h, k):
    return np.asarray(h[k][:]) if k in h else None
def neq(a, b):
    if a is None or b is None: return None
    if a.shape != b.shape: return -1
    return int(np.count_nonzero(a != b))
def reldiff(a, b, floor_rel=1e-30):
    a = np.asarray(a, np.float64); b = np.asarray(b, np.float64)
    den = np.maximum(np.abs(b), floor_rel * np.max(np.abs(b)) if b.size else 1.0)
    return np.abs(a - b) / den

print("== 0. ダンプの属性")
F = {}
for c, (p32, p64) in DUMPS.items():
    try:
        F[c] = (h5py.File(p32 / "stage2.h5", "r"), h5py.File(p64 / "stage2.h5", "r"))
    except Exception as e:
        print(f"  {c}: ダンプが読めない ({e})"); continue
    a32, a64 = F[c][0].attrs, F[c][1].attrs
    print(f"  {c}: float {a32['flow_float_bytes']} B / FP64 {a64['flow_float_bytes']} B, geom64 {a32['has_geom64']}/{a64['has_geom64']}, "
          f"nCells {a32['nCells']}, gradLSQ {a32['gradLSQ']}, lineImplicit {a32['lineImplicit']}, axisym {a32['isAxisymmetric']}, "
          f"lines device==new {a32.get('lines_device_equal_new')}/{a64.get('lines_device_equal_new')}, "
          f"LSQ production==rebuild mismatch {a32.get('lsq_production_vs_rebuild_mismatch')}/{a64.get('lsq_production_vs_rebuild_mismatch')}")
    print(f"     CSR 一致 float vs FP64: cell_planes {neq(ds(F[c][0], '/mesh/cell_planes'), ds(F[c][1], '/mesh/cell_planes'))} 件不一致, "
          f"plane_cells {neq(ds(F[c][0], '/mesh/plane_cells'), ds(F[c][1], '/mesh/plane_cells'))} 件不一致")

print("== 1. 接続の一致 (新しい float 対 FP64。() は旧 float 対 FP64)")
for c in ("case/45", "case/56"):
    if c not in F: verdict(f"1.lines.{c}", False, "ダンプなし"); continue
    f32, f64 = F[c]
    mis = sum(neq(ds(f32, f"/lines/new/{k}"), ds(f64, f"/lines/new/{k}")) or 0 for k in ("offsets", "cells"))
    old = {k: neq(ds(f32, f"/lines/legacy/{k}"), ds(f64, f"/lines/new/{k}")) for k in ("offsets", "cells", "prev", "next")}
    dev = (f32.attrs.get("lines_device_equal_new"), f64.attrs.get("lines_device_equal_new"))
    verdict(f"1.lines.{c}", mis == 0, f"offsets・cells の不一致 {mis} (旧 float: {old})、デバイス == new {dev}")
if "case/39" in F:
    f32, f64 = F["case/39"]
    ids = [k for k in f64["/periodic"].keys() if k != "periodicRoot"] if "/periodic" in f64 else []
    mis = sum(neq(ds(f32, f"/periodic/{i}/partnerPlnID"), ds(f64, f"/periodic/{i}/partnerPlnID")) or 0 for i in ids)
    old = sum(neq(ds(f32, f"/periodic/{i}/partnerPlnID_legacy"), ds(f64, f"/periodic/{i}/partnerPlnID")) or 0 for i in ids)
    root = neq(ds(f32, "/periodic/periodicRoot"), ds(f64, "/periodic/periodicRoot"))
    verdict("1.periodic.case/39", mis == 0 and root == 0 and len(ids) > 0, f"周期 bcond {ids}、相手の不一致 {mis} (旧 float {old})、periodicRoot の不一致 {root}")
    sm = sum(neq(ds(f32, f"/lsq/new/{k}"), ds(f64, f"/lsq/new/{k}")) or 0 for k in ("seam_group", "seam_class"))
    so = sum(neq(ds(f32, f"/lsq/legacy/{k}"), ds(f64, f"/lsq/new/{k}")) or 0 for k in ("seam_group", "seam_class"))
    ng = len(np.unique(ds(f64, "/lsq/new/seam_group")[ds(f64, "/lsq/new/seam_group") >= 0])) if "/lsq/new/seam_group" in f64 else 0
    verdict("1.seam.case/39", sm == 0 and ng > 0, f"継ぎ目の group {ng}、同値類の不一致 {sm} (旧 float {so})")
for c in ("case/48", "case/39", "case/45", "case/56"):
    if c not in F or "/wall" not in F[c][1]: continue
    f32, f64 = F[c]
    tot = mis = old = 0
    for i in f64["/wall"].keys():
        a, b, l = ds(f32, f"/wall/{i}/irep"), ds(f64, f"/wall/{i}/irep"), ds(f32, f"/wall/{i}/irep_legacy")
        tot += b.size; mis += neq(a, b) or 0; old += neq(l, b) or 0
    key = "1.wall_irep.case/48" if c == "case/48" else f"1.wall_irep.{c} (記録)"
    verdict(key, mis == 0 and tot > 0, f"境界面 {tot}、irep の不一致 {mis} (旧 float {old})")

print("== 2. 係数の精度 (新しい float 対 FP64。() は旧 float)")
for c in DUMPS:
    if c not in F: continue
    f32, f64 = F[c]
    if "/lsq/new/cInt" in f64:
        c64 = ds(f64, "/lsq/new/cInt").astype(np.float64)
        rn = reldiff(ds(f32, "/lsq/new/cInt"), c64); ro = reldiff(ds(f32, "/lsq/legacy/cInt"), c64)
        p = np.percentile(rn, 99.9)
        verdict(f"2.cInt.{c}", p <= 1e-5, f"相対差 99.9 % 点 {p:.3g}・最大 {rn.max():.3g} (旧 99.9 % {np.percentile(ro, 99.9):.3g}・最大 {ro.max():.3g}・>1e-3 {int(np.sum(ro > 1e-3))}/{ro.size})")
    if "/wall" in f64:
        worst = worst_o = 0.0; n = 0
        for i in f64["/wall"].keys():
            ok = (ds(f64, f"/wall/{i}/irep") >= 0) & (ds(f32, f"/wall/{i}/irep") == ds(f64, f"/wall/{i}/irep"))
            y64 = np.maximum(ds(f64, f"/wall/{i}/dn").astype(np.float64), 1e-12)[ok]
            y32 = np.maximum(ds(f32, f"/wall/{i}/dn").astype(np.float64), 1e-12)[ok]
            okl = ok & (ds(f32, f"/wall/{i}/irep_legacy") == ds(f64, f"/wall/{i}/irep"))
            yl = np.maximum(ds(f32, f"/wall/{i}/dn_legacy").astype(np.float64), 1e-12)[okl]
            yl64 = np.maximum(ds(f64, f"/wall/{i}/dn").astype(np.float64), 1e-12)[okl]
            if y64.size: worst = max(worst, float(np.max(np.abs(y32 - y64) / y64)))
            if yl.size: worst_o = max(worst_o, float(np.max(np.abs(yl - yl64) / yl64)))
            n += int(ok.sum())
        verdict(f"2.wall_y.{c}", worst <= 1e-6, f"y の相対差の最大 {worst:.3g} (面 {n}) (旧 float {worst_o:.3g})")
    if "/d1d2" in f64:
        worst = worst_o = 0.0; n = 0; jm = 0
        for i in f64["/d1d2"].keys():
            for k in ("d1", "d2"):
                a, b, l = ds(f32, f"/d1d2/{i}/{k}"), ds(f64, f"/d1d2/{i}/{k}"), ds(f32, f"/d1d2/{i}/{k}_legacy")
                m = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
                if m.any(): worst = max(worst, float(np.max(np.abs(a[m] - b[m]) / np.abs(b[m]))))
                ml = np.isfinite(l) & np.isfinite(b) & (np.abs(b) > 0)
                if ml.any(): worst_o = max(worst_o, float(np.max(np.abs(l[ml] - b[ml]) / np.abs(b[ml]))))
                n += int(m.sum())
            jm += (neq(ds(f32, f"/d1d2/{i}/jdof"), ds(f64, f"/d1d2/{i}/jdof")) or 0) + (neq(ds(f32, f"/d1d2/{i}/jdof2"), ds(f64, f"/d1d2/{i}/jdof2")) or 0)
        verdict(f"2.d1d2.{c}", worst <= 1e-6 and jm == 0, f"d1・d2 の相対差の最大 {worst:.3g} (値 {n})、jdof の不一致 {jm} (旧 float {worst_o:.3g})")
    if "/delta_les" in f64:
        nc = int(f64.attrs["nCells"]); b = ds(f64, "/delta_les")[:nc].astype(np.float64); m = b > 0
        rn = np.abs(ds(f32, "/delta_les")[:nc][m] - b[m]) / b[m]; ro = np.abs(ds(f32, "/delta_les_legacy")[:nc][m] - b[m]) / b[m]
        print(f"  (記録) delta_les.{c}: 相対差の最大 {rn.max():.3g} (旧 float {ro.max():.3g})")
    if "/axisym/A_closure_y" in f64:
        nc = int(f64.attrs["nCells"])
        for ax in ("x", "y"):
            a, b = ds(f32, f"/axisym/A_closure_{ax}")[:nc].astype(np.float64), ds(f64, f"/axisym/A_closure_{ax}")[:nc].astype(np.float64)
            s = ds(f64, f"/axisym/closure_abs64_{ax}")[:nc]; m = s > 0
            print(f"  (記録) A_closure_{ax}.{c}: |float − FP64| / Σ|S| の最大 {np.max(np.abs(a[m] - b[m]) / s[m]):.3g}・RMS {np.sqrt(np.mean((np.abs(a[m] - b[m]) / s[m]) ** 2)):.3g}")

print("== 3. closure の非劣化 (一様な静止場、P0 = 1e5 Pa)")
def vals(r, n, keys=None):
    with h5py.File(r / f"res_{n}.h5", "r") as h:
        return {k: np.asarray(h["VALUE/" + k][:], dtype=np.float64) for k in h["VALUE"].keys() if keys is None or k in keys}
P0 = 1.0e5
S = {k: C45 / f"run_{k}" for k in ("0399_fg2_still_old32", "0400_fg2_still_new32", "0401_fg2_still_old64", "0402_fg2_still_new64", "0403_fg2_still_old64b", "0404_fg2_still_old32b")}
try:
    Ap = {32: ds(F["case/45"][0], "/axisym/A_planar").astype(np.float64), 64: ds(F["case/45"][1], "/axisym/A_planar").astype(np.float64)}
    st = {}
    for k, r in S.items():
        prec = 32 if "32" in k.split("_")[-1] else 64
        v1 = vals(r, 1, ("res_roUy", "P", "Ux", "Uy", "Uz")); v10 = vals(r, 10, ("Ux", "Uy", "Uz", "P", "ro", "T"))
        n = v1["res_roUy"].size; q = v1["res_roUy"] / (P0 * Ap[prec][:n])
        bad = sum(int(np.count_nonzero(~np.isfinite(x))) for x in list(v1.values()) + list(v10.values()))
        um = float(np.max(np.sqrt(v10["Ux"] ** 2 + v10["Uy"] ** 2 + v10["Uz"] ** 2)))
        st[k] = dict(q=q, qmax=float(np.max(np.abs(q))), qrms=float(np.sqrt(np.mean(q ** 2))), umax=um, nonfinite=bad,
                     Pdev=float(np.max(np.abs(v1["P"] - P0)) / P0))
        j = np.arange(n) % 121; im = int(np.argmax(np.abs(q))); inn = (j >= 1) & (j <= 119)
        u10 = np.sqrt(v10["Ux"] ** 2 + v10["Uy"] ** 2 + v10["Uz"] ** 2); iu = int(np.argmax(u10))
        print(f"  {k}: |res_roUy|/(P0 A_planar) 最大 {st[k]['qmax']:.3e} (節点 {im}: i {im // 121}・j {im % 121})・RMS {st[k]['qrms']:.3e}"
              f" [軸と壁を除く 1≤j≤119: 最大 {np.max(np.abs(q[inn])):.3e}・RMS {np.sqrt(np.mean(q[inn] ** 2)):.3e}]、"
              f"10 step 後の最大 |u| {um:.3e} m/s (i {iu // 121}・j {iu % 121})、step 1 の |P − P0|/P0 最大 {st[k]['Pdev']:.2e}、非有限 {bad}")
    o, nw = st["0399_fg2_still_old32"], st["0400_fg2_still_new32"]
    verdict("3.closure.float.res", nw["qmax"] <= 1.1 * o["qmax"] and nw["qrms"] <= 1.1 * o["qrms"],
            f"新 / 旧: 最大 {nw['qmax'] / o['qmax']:.3f}・RMS {nw['qrms'] / o['qrms']:.3f} (旧の再実行 / 旧: 最大 {st['0404_fg2_still_old32b']['qmax'] / o['qmax']:.3f}・RMS {st['0404_fg2_still_old32b']['qrms'] / o['qrms']:.3f})")
    verdict("3.closure.float.u10", nw["umax"] <= 1.1 * o["umax"],
            f"新 / 旧: {nw['umax'] / o['umax']:.3f} (旧の再実行 / 旧 {st['0404_fg2_still_old32b']['umax'] / o['umax']:.3f})")
    a, b, c = st["0401_fg2_still_old64"]["q"], st["0402_fg2_still_new64"]["q"], st["0403_fg2_still_old64b"]["q"]
    dn, dr = float(np.max(np.abs(b - a))), float(np.max(np.abs(c - a)))
    verdict("3.closure.fp64", dn <= 10 * dr if dr > 0 else dn == 0, f"新−旧の最大 {dn:.3e}、再実行の差 {dr:.3e} (比 {dn / dr if dr > 0 else float('nan'):.2f})")
except Exception as e:
    verdict("3.closure", False, f"集計できない ({e})")

print("== 4・5. 既定の経路の同一性 (20 step 後の場: 最大 |差| / 最大 |値|)")
def hist(r):
    rows = [x for x in csv.DictReader(open(r / "residual_history.csv")) if x["phase"] == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
    return np.array([[float(x[c]) for c in cols] for x in rows])
for key, base, (o, nw, o2, s1) in (("4.fp64.case/45", C45, ("run_0390_fg_id64_old", "run_0398_fg2_id64_new", "run_0392_fg_id64_old2", "run_0391_fg_id64_new")),
                                   ("5.float.case/56", C56, ("run_0079_fg_id32_old", "run_0084_fg2_id32_new", "run_0081_fg_id32_old2", "run_0080_fg_id32_new"))):
    try:
        A, B, A2, S1 = (vals(base / x, 20) for x in (o, nw, o2, s1))
        worst = 0.0; lines = []
        for k in ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T"):
            if k not in A: continue
            s = np.max(np.abs(A[k])); dn = np.max(np.abs(B[k] - A[k])) / s; dr = np.max(np.abs(A2[k] - A[k])) / s; d1 = np.max(np.abs(B[k] - S1[k])) / s
            ratio = dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf); worst = max(worst, ratio)
            lines.append(f"{k} {dn:.2e}/{dr:.2e} (段①比 {d1:.2e})")
        nonf = sum(int(np.count_nonzero(~np.isfinite(B[k]))) for k in B)
        ha, hb, ha2 = hist(base / o), hist(base / nw), hist(base / o2)
        rn = np.max(np.abs(hb - ha) / np.maximum(np.abs(ha), 1e-300)); rr = np.max(np.abs(ha2 - ha) / np.maximum(np.abs(ha), 1e-300))
        verdict(key, worst <= 10 and nonf == 0, f"新−旧 / 再実行の差の最大比 {worst:.2f}、非有限 {nonf}; " + "; ".join(lines) + f"; 残差履歴 新−旧 {rn:.2e} / 再実行 {rr:.2e}")
    except Exception as e:
        verdict(key, False, f"集計できない ({e})")

json.dump(V, open(C45 / "fg2_verdict.json", "w"), ensure_ascii=False, indent=1)
print("== 総合:", "ALL PASS" if all(v["pass"] for k, v in V.items() if "(記録)" not in k) else "FAIL あり",
      [k for k, v in V.items() if not v["pass"]])
