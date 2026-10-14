"""plan architecture-float-state-double-geometry §6.5 の 4a・4b・4c (2026-10-10、codex diagnose の採用) — 段 ② の判定のやり直し。
fg2_an.py (最初の集計、記録として残す) の欠陥を直す: 欠損・形の不一致・非有限・本番の配列との不一致は判定不能 (= FAIL)。
4a は LSQ の係数の改めた基準 (ベクトルの相対差・有限・ゼロベクトル・線形の再現の行列)、4b は closure を使う経路。
AWS の case/45 で回す。結果を標準出力と fg2_rejudge.json に書く。"""
import csv, json, h5py, numpy as np
from pathlib import Path
H = Path("/home/ubuntu/forge-wallfit/case")
C45, C56, C39, C48 = H / "45.isobutane_m6_d155", H / "56.gap_tp1187", H / "39.periodic_hills", H / "48.flat_plate_cooled_m4"
DUMPS = {"case/45": (C45 / "run_0396_fg2_dump32", C45 / "run_0397_fg2_dump64", C45 / "run_0402_fg2_still_new64/res_0.h5"),
         "case/56": (C56 / "run_0082_fg2_dump32", C56 / "run_0083_fg2_dump64", C56 / "run_0082_fg2_dump32/mesh.h5"),
         "case/39": (C39 / "run_0956_fg2_dump32", C39 / "run_0957_fg2_dump64", C39 / "run_0956_fg2_dump32/hill_xc_80x50x30.h5"),
         "case/48": (C48 / "run_0046_fg2_dump32", C48 / "run_0047_fg2_dump64", C48 / "run_0046_fg2_dump32/mesh.h5")}
V = {}
class Undecidable(Exception): pass
def verdict(key, ok, detail):
    V[key] = {"verdict": "PASS" if ok else "FAIL", "detail": detail}
    print(f"  [{'PASS' if ok else 'FAIL'}] {key}: {detail}")
def undecidable(key, why):
    V[key] = {"verdict": "UNDECIDABLE", "detail": why}
    print(f"  [判定不能] {key}: {why}")
def ds(h, k):
    if k not in h: raise Undecidable(f"{h.filename}: {k} が無い")
    a = np.asarray(h[k][:])
    if a.dtype.kind == "f" and not np.all(np.isfinite(a)): raise Undecidable(f"{h.filename}: {k} に非有限 {int(np.count_nonzero(~np.isfinite(a)))}")
    return a
def neq(a, b, what):
    if a.shape != b.shape: raise Undecidable(f"{what}: 形が違う {a.shape} / {b.shape}")
    return int(np.count_nonzero(a != b))
def attr(h, k):
    if k not in h.attrs: raise Undecidable(f"{h.filename}: 属性 {k} が無い")
    return h.attrs[k]
def judge(key, fn):
    try: fn()
    except Undecidable as e: undecidable(key, str(e))

F = {c: (h5py.File(a / "stage2.h5", "r"), h5py.File(b / "stage2.h5", "r")) for c, (a, b, _) in DUMPS.items()}

print("== 4c 本番の配列との一致 (ダンプの属性)")
for c, (f32, f64) in F.items():
    def chk(c=c, f32=f32, f64=f64):
        bad = []
        for f in (f32, f64):
            if int(attr(f, "lineImplicit")) == 1 and int(attr(f, "lines_device_equal_new")) != 1: bad.append(f"{f.filename}: lines_device_equal_new ≠ 1")
            if int(attr(f, "gradLSQ")) == 2 and int(attr(f, "lsq_production_vs_rebuild_mismatch")) != 0: bad.append(f"{f.filename}: LSQ の本番と作り直しが不一致")
            if int(attr(f, "has_geom64")) != 1: bad.append(f"{f.filename}: has_geom64 ≠ 1")
        verdict(f"4c.production.{c}", not bad, "; ".join(bad) or "ライン (lineImplicit のとき)・LSQ の本番の配列 = ダンプの新しい作り方、double の写しあり")
    judge(f"4c.production.{c}", chk)

print("== 1 接続の一致 (新しい float 対 FP64。欠損は判定不能)")
for c in ("case/45", "case/56"):
    def f1(c=c):
        f32, f64 = F[c]
        mis = sum(neq(ds(f32, f"/lines/new/{k}"), ds(f64, f"/lines/new/{k}"), f"{c} lines {k}") for k in ("offsets", "cells", "prev", "next"))
        old = sum(neq(ds(f32, f"/lines/legacy/{k}"), ds(f64, f"/lines/new/{k}"), f"{c} lines legacy {k}") for k in ("offsets", "cells", "prev", "next"))
        verdict(f"1.lines.{c}", mis == 0, f"offsets・cells・prev・next の不一致 {mis} (旧 float {old})")
    judge(f"1.lines.{c}", f1)
def f1p():
    f32, f64 = F["case/39"]
    ids = sorted(k for k in f64["/periodic"].keys() if k != "periodicRoot")
    if not ids: raise Undecidable("周期 bcond が無い")
    mis = sum(neq(ds(f32, f"/periodic/{i}/{k}"), ds(f64, f"/periodic/{i}/{k}"), f"periodic {i} {k}") for i in ids for k in ("partnerPlnID", "partnerCellID"))
    old = sum(neq(ds(f32, f"/periodic/{i}/partnerPlnID_legacy"), ds(f64, f"/periodic/{i}/partnerPlnID"), f"periodic {i} legacy") for i in ids)
    root = neq(ds(f32, "/periodic/periodicRoot"), ds(f64, "/periodic/periodicRoot"), "periodicRoot")
    verdict("1.periodic.case/39", mis == 0 and root == 0, f"bcond {ids}、相手 (面・CV) の不一致 {mis} (旧 float {old})、periodicRoot の不一致 {root}")
    g64 = ds(f64, "/lsq/new/seam_group"); ng = len(np.unique(g64[g64 >= 0]))
    if ng == 0: raise Undecidable("継ぎ目の group が 0")
    sm = sum(neq(ds(f32, f"/lsq/new/{k}"), ds(f64, f"/lsq/new/{k}"), f"seam {k}") for k in ("seam_group", "seam_class"))
    so = sum(neq(ds(f32, f"/lsq/legacy/{k}"), ds(f64, f"/lsq/new/{k}"), f"seam legacy {k}") for k in ("seam_group", "seam_class"))
    verdict("1.seam.case/39", sm == 0, f"group {ng}、同値類の不一致 {sm} (旧 float {so})")
judge("1.periodic.case/39", f1p)
for c in ("case/48", "case/39", "case/45", "case/56"):
    def f1w(c=c):
        f32, f64 = F[c]
        ids = sorted(f64["/wall"].keys()) if "/wall" in f64 else []
        if not ids: raise Undecidable("壁 bcond が無い")
        tot = sum(ds(f64, f"/wall/{i}/irep").size for i in ids)
        mis = sum(neq(ds(f32, f"/wall/{i}/irep"), ds(f64, f"/wall/{i}/irep"), f"irep {i}") for i in ids)
        old = sum(neq(ds(f32, f"/wall/{i}/irep_legacy"), ds(f64, f"/wall/{i}/irep"), f"irep legacy {i}") for i in ids)
        nocand = sum(int(np.count_nonzero(ds(f64, f"/wall/{i}/irep") < 0)) for i in ids)
        verdict(f"1.wall_irep.{c}" + ("" if c == "case/48" else " (記録)"), mis == 0, f"境界面 {tot}、irep の不一致 {mis} (旧 float {old})、候補なし {nocand}")
    judge(f"1.wall_irep.{c}", f1w)

print("== 2・4a 係数の精度")
for c, (_, _, meshf) in DUMPS.items():
    def f4a(c=c, meshf=meshf):
        f32, f64 = F[c]
        nC = int(attr(f64, "nCells")); nN = int(attr(f64, "nNormalPlanes"))
        c32 = ds(f32, "/lsq/new/cInt").astype(np.float64).reshape(-1, 3); c64 = ds(f64, "/lsq/new/cInt").astype(np.float64).reshape(-1, 3)
        co = ds(f32, "/lsq/legacy/cInt").astype(np.float64).reshape(-1, 3)
        if c32.shape != c64.shape: raise Undecidable(f"cInt の形が違う {c32.shape} / {c64.shape}")
        cpi = ds(f64, "/mesh/cell_planes_index"); cp = ds(f64, "/mesh/cell_planes"); pc = ds(f64, "/mesh/plane_cells").reshape(-1, 2)
        if neq(cpi, ds(f32, "/mesh/cell_planes_index"), "cpi") or neq(cp, ds(f32, "/mesh/cell_planes"), "cp"): raise Undecidable("CSR が float と FP64 で違う")
        nInc = int(cpi[-1])
        if nInc != c64.shape[0]: raise Undecidable(f"incidence の数が合わない {nInc} / {c64.shape[0]}")
        owner = np.repeat(np.arange(nC), np.diff(cpi[:nC + 1]))
        ip = cp.astype(np.int64); internal = ip < nN
        # 元の判定 (記録): 係数ごとの相対差の 99.9 % 点
        mx = np.max(np.abs(c64)); r_comp = np.abs(c32 - c64) / np.maximum(np.abs(c64), 1e-30 * mx)
        n64 = np.linalg.norm(c64, axis=1); zero64 = n64 == 0; zero32 = np.linalg.norm(c32, axis=1) == 0
        nz = ~zero64
        rv = np.linalg.norm(c32 - c64, axis=1)[nz] / n64[nz]; rvo = np.linalg.norm(co - c64, axis=1)[nz] / n64[nz]
        zero_ok = bool(np.all(zero32[zero64])); bnd_ok = bool(np.all(zero64[~internal]) and np.all(zero32[~internal]))
        # 線形の再現の行列 M_i = Σ c ⊗ d (d は同じ倍精度の変位を両方に使う)
        with h5py.File(meshf, "r") as m: X = np.asarray(m["MESH/COORD"][:], np.float64).reshape(-1, 3)
        if X.shape[0] < nC: raise Undecidable(f"座標の数 {X.shape[0]} < nCells {nC}")
        a, b = pc[ip.clip(0, len(pc) - 1), 0], pc[ip.clip(0, len(pc) - 1), 1]
        jc = np.where(a == owner, b, a); d = np.where(internal[:, None], X[np.clip(jc, 0, nC - 1)] - X[owner], 0.0)
        M32 = np.zeros((nC, 3, 3)); M64 = np.zeros((nC, 3, 3)); Mo = np.zeros((nC, 3, 3))
        for k in range(3):
            for l in range(3):
                M32[:, k, l] = np.bincount(owner, weights=c32[:, k] * d[:, l], minlength=nC)
                M64[:, k, l] = np.bincount(owner, weights=c64[:, k] * d[:, l], minlength=nC)
                Mo[:, k, l] = np.bincount(owner, weights=co[:, k] * d[:, l], minlength=nC)
        f64n = np.linalg.norm(M64.reshape(nC, -1), axis=1); okn = f64n > 0
        rM = np.linalg.norm((M32 - M64).reshape(nC, -1), axis=1)[okn] / f64n[okn]
        rMo = np.linalg.norm((Mo - M64).reshape(nC, -1), axis=1)[okn] / f64n[okn]
        tr = np.trace(M64, axis1=1, axis2=2)[okn]
        ok = rv.max() <= 1e-6 and zero_ok and bnd_ok and rM.max() <= 1e-6
        print(f"  (記録) 2.cInt.{c} 元の判定 (係数ごとの 99.9 % 点 ≤ 1e-5): {np.percentile(r_comp, 99.9):.3g} → {'PASS' if np.percentile(r_comp, 99.9) <= 1e-5 else 'FAIL'}")
        verdict(f"4a.cInt.{c}", ok,
                f"ベクトルの相対差の最大 {rv.max():.3g} (旧 float {rvo.max():.3g})、ゼロベクトル {int(zero64.sum())} 本で float も 0: {zero_ok}、"
                f"境界の incidence {int((~internal).sum())} 本が両方 0: {bnd_ok}、線形の再現の行列の差の最大 {rM.max():.3g} (旧 float {rMo.max():.3g}; "
                f"tr M64 の範囲 {tr.min():.4f}〜{tr.max():.4f}、節点 {int(okn.sum())})")
    judge(f"4a.cInt.{c}", f4a)
    def f2y(c=c):
        f32, f64 = F[c]
        ids = sorted(f64["/wall"].keys()) if "/wall" in f64 else []
        if not ids: raise Undecidable("壁 bcond が無い")
        worst = worst_o = 0.0; n = 0
        for i in ids:
            irep64 = ds(f64, f"/wall/{i}/irep")
            ok = (irep64 >= 0) & (ds(f32, f"/wall/{i}/irep") == irep64)
            y64 = np.maximum(ds(f64, f"/wall/{i}/dn").astype(np.float64), 1e-12); y32 = np.maximum(ds(f32, f"/wall/{i}/dn").astype(np.float64), 1e-12)
            yl = np.maximum(ds(f32, f"/wall/{i}/dn_legacy").astype(np.float64), 1e-12); okl = ok & (ds(f32, f"/wall/{i}/irep_legacy") == irep64)
            if ok.any(): worst = max(worst, float(np.max(np.abs(y32[ok] - y64[ok]) / y64[ok])))
            if okl.any(): worst_o = max(worst_o, float(np.max(np.abs(yl[okl] - y64[okl]) / y64[okl])))
            n += int(ok.sum())
        if n == 0: raise Undecidable("比べられる面が 0")
        verdict(f"2.wall_y.{c}", worst <= 1e-6, f"y の相対差の最大 {worst:.3g} (面 {n}) (旧 float {worst_o:.3g})")
    judge(f"2.wall_y.{c}", f2y)
    def f2d(c=c):
        f32, f64 = F[c]
        ids = sorted(f64["/d1d2"].keys()) if "/d1d2" in f64 else []
        if not ids: raise Undecidable("d1d2 が無い")
        worst = 0.0; n = 0; jm = 0; okm = 0
        for i in ids:
            okk = ds(f64, f"/d1d2/{i}/ok").astype(bool)
            okm += neq(ds(f32, f"/d1d2/{i}/ok"), ds(f64, f"/d1d2/{i}/ok"), f"d1d2 ok {i}")
            for k in ("jdof", "jdof2"): jm += neq(ds(f32, f"/d1d2/{i}/{k}"), ds(f64, f"/d1d2/{i}/{k}"), f"d1d2 {k} {i}")
            for k in ("d1", "d2"):
                a = np.asarray(f32[f"/d1d2/{i}/{k}"][:], np.float64)[okk]; b = np.asarray(f64[f"/d1d2/{i}/{k}"][:], np.float64)[okk]
                if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))): raise Undecidable(f"d1d2 {i} {k}: ok の面に非有限")
                m = np.abs(b) > 0; n += int(m.sum())
                if m.any(): worst = max(worst, float(np.max(np.abs(a[m] - b[m]) / np.abs(b[m]))))
        verdict(f"2.d1d2.{c}", worst <= 1e-6 and jm == 0 and okm == 0, f"ok の面の d1・d2 の相対差の最大 {worst:.3g} (値 {n})、jdof の不一致 {jm}、ok の不一致 {okm}")
    judge(f"2.d1d2.{c}", f2d)

print("== 4b closure を使う経路 (hoopAreaFromClosure 1 のダンプ)")
def f4b():
    try:
        g32 = h5py.File(C45 / "run_0415_fg2_dumphc32/stage2.h5", "r"); g64 = h5py.File(C45 / "run_0416_fg2_dumphc64/stage2.h5", "r")
    except OSError as e:
        raise Undecidable(f"closure を使うダンプが無い ({e})")
    msg = []; ok = True
    for g, tag in ((g32, "float"), (g64, "FP64")):
        if int(attr(g, "closure_active")) != 1: raise Undecidable(f"{tag}: closure_active ≠ 1")
        nc = int(attr(g, "nCells")); dt = np.float32 if int(attr(g, "flow_float_bytes")) == 4 else np.float64
        for ax in ("x", "y"):
            A = ds(g, f"/axisym/A_closure_{ax}")[:nc]; S = ds(g, f"/axisym/closure_sum64_{ax}")[:nc].astype(dt)
            mis = int(np.count_nonzero(A != S)); ok &= mis == 0
            msg.append(f"{tag} {ax}: デバイスの A_closure と Σ±S (double の和を 1 回丸めたもの) の不一致 {mis}/{nc}")
    nc = int(attr(g64, "nCells"))
    for ax in ("x", "y"):
        a, b = ds(g32, f"/axisym/A_closure_{ax}")[:nc].astype(np.float64), ds(g64, f"/axisym/A_closure_{ax}")[:nc].astype(np.float64)
        s = ds(g64, f"/axisym/closure_abs64_{ax}")[:nc]; m = s > 0
        msg.append(f"(記録) float − FP64 の A_closure_{ax} / Σ|S| の最大 {np.max(np.abs(a[m] - b[m]) / s[m]):.3g}")
        raw = ds(g64, f"/axisym/closure_raw64_{ax}")[:nc]
        msg.append(f"(記録) FP64 の A_closure_{ax} と生の double の幾何からの closure の差 / Σ|S| の最大 {np.max(np.abs(b[m] - raw[m]) / s[m]):.3g}")
    verdict("4b.closure_sum", ok, "; ".join(msg))
judge("4b.closure_sum", f4b)
def vals(r, n, keys):
    with h5py.File(r / f"res_{n}.h5", "r") as h:
        out = {k: np.asarray(h["VALUE/" + k][:], dtype=np.float64) for k in keys}
    for k, v in out.items():
        if not np.all(np.isfinite(v)): raise Undecidable(f"{r.name} res_{n} の {k} に非有限")
    return out
def f4b_fp64():
    A = ds(F["case/45"][1], "/axisym/A_planar").astype(np.float64)
    q = {}
    for r in ("run_0408_fg2_stillhc_new64", "run_0417_fg2_stillhc_old64", "run_0418_fg2_stillhc_old64b"):
        p = C45 / r
        if not (p / "res_1.h5").exists(): raise Undecidable(f"{r} の res_1 が無い")
        v = vals(p, 1, ("res_roUy",)); vals(p, 10, ("Ux", "Uy", "P", "T", "ro"))
        q[r] = v["res_roUy"] / (1.0e5 * A[:v["res_roUy"].size])
    dn = float(np.max(np.abs(q["run_0408_fg2_stillhc_new64"] - q["run_0417_fg2_stillhc_old64"])))
    dr = float(np.max(np.abs(q["run_0418_fg2_stillhc_old64b"] - q["run_0417_fg2_stillhc_old64"])))
    verdict("4b.closure_fp64_identity", dn <= 10 * dr if dr > 0 else dn == 0, f"静止場 (closure あり) の FP64 の新−旧の最大 {dn:.3e}、再実行の差 {dr:.3e}")
judge("4b.closure_fp64_identity", f4b_fp64)

print("== 3・3' (非有限を合否に入れる)")
def f3(prefix, runs, key):
    def fn():
        A = {32: ds(F["case/45"][0], "/axisym/A_planar").astype(np.float64), 64: ds(F["case/45"][1], "/axisym/A_planar").astype(np.float64)}
        st = {}
        for r in runs:
            prec = 64 if r.endswith("64") or r.endswith("64b") else 32
            v1 = vals(C45 / r, 1, ("res_roUy",)); v10 = vals(C45 / r, 10, ("Ux", "Uy", "Uz", "P", "T", "ro"))
            q = v1["res_roUy"] / (1.0e5 * A[prec][:v1["res_roUy"].size])
            st[r] = (float(np.max(np.abs(q))), float(np.sqrt(np.mean(q ** 2))), float(np.max(np.sqrt(v10["Ux"] ** 2 + v10["Uy"] ** 2 + v10["Uz"] ** 2))))
        o, n, ob = (st[x] for x in runs[:3])
        verdict(key, n[0] <= 1.1 * o[0] and n[1] <= 1.1 * o[1] and n[2] <= 1.1 * o[2],
                f"新 / 旧: 最大 {n[0] / o[0]:.3f}・RMS {n[1] / o[1]:.3f}・|u| {n[2] / o[2]:.3f} (旧の再実行 / 旧 {ob[0] / o[0]:.3f}・{ob[1] / o[1]:.3f}・{ob[2] / o[2]:.3f})、非有限なし")
    judge(key, fn)
f3("3", ["run_0399_fg2_still_old32", "run_0400_fg2_still_new32", "run_0404_fg2_still_old32b"], "3.closure.float (既定の経路)")
f3("3p", ["run_0405_fg2_stillhc_old32", "run_0406_fg2_stillhc_new32", "run_0407_fg2_stillhc_old32b"], "3'.closure.float (closure あり)")

print("== 4・5 既定の経路の同一性 (非有限を合否に入れる)")
def hist(r):
    rows = [x for x in csv.DictReader(open(r / "residual_history.csv")) if x["phase"] == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
    a = np.array([[float(x[c]) for c in cols] for x in rows])
    if not np.all(np.isfinite(a)): raise Undecidable(f"{r.name} の残差履歴に非有限")
    return a
for key, base, (o, nw, o2) in (("4.fp64.case/45", C45, ("run_0390_fg_id64_old", "run_0398_fg2_id64_new", "run_0392_fg_id64_old2")),
                               ("5.float.case/56", C56, ("run_0079_fg_id32_old", "run_0084_fg2_id32_new", "run_0081_fg_id32_old2"))):
    def f45(key=key, base=base, o=o, nw=nw, o2=o2):
        ks = ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T")
        A, B, A2 = (vals(base / x, 20, ks) for x in (o, nw, o2))
        worst = 0.0
        for k in ks:
            s = np.max(np.abs(A[k])); dn = np.max(np.abs(B[k] - A[k])) / s; dr = np.max(np.abs(A2[k] - A[k])) / s
            worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf))
        ha, hb, ha2 = hist(base / o), hist(base / nw), hist(base / o2)
        verdict(key, worst <= 10, f"新−旧 / 再実行の差の最大比 {worst:.2f}、場・残差履歴に非有限なし")
    judge(key, f45)

json.dump(V, open(C45 / "fg2_rejudge.json", "w"), ensure_ascii=False, indent=1)
bad = [k for k, v in V.items() if v["verdict"] != "PASS" and "(記録)" not in k]
print("== 総合:", "ALL PASS" if not bad else f"PASS でないもの {bad}")
