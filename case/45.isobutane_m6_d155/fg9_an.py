"""plan axisymmetric-freestream-hoop-gauge §4.6 (事前登録) の 2〜6 の判定。fg9.sh の後に AWS の case/45 で回す。結果は fg9_verdict.json。
欠損・非有限は判定不能 (= PASS でない)。"""
import json, re, h5py, numpy as np
from pathlib import Path
C45 = Path("/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"); C48 = C45.parent / "48.flat_plate_cooled_m4"; C56 = C45.parent / "56.gap_tp1187"
W = C45 / "_conv9"; NJ = 121; RT = 0.076807; EPS = np.finfo(np.float64).eps; P0 = 1.0e5
V = {}
def v(key, ok, det): V[key] = {"verdict": "PASS" if ok else "FAIL", "detail": det}; print(f"  [{'PASS' if ok else 'FAIL'}] {key}: {det}")
def und(key, why): V[key] = {"verdict": "UNDECIDABLE", "detail": why}; print(f"  [判定不能] {key}: {why}")
log = (C45 / "fg9.log").read_text()   # fg9.sh と fg9b.sh の両方の記録 (正規表現は最後の一致を使う)
print("== 2 変換器")
last = lambda pat: (list(re.finditer(pat, log)) or [None])[-1]
m = last(r"interior faces: rebuilt W vs rSurfVect max \|ΔW\|/\|W\| ([0-9.eE+-]+)"); mb = last(r"boundary half faces \((\d+), unmapped (\d+)\): rebuilt W vs rSurfVect max \|ΔW\|/\|W\| ([0-9.eE+-]+)")
if m and mb: v("2a.rSurfVect", float(m.group(1)) <= 1e-12 and float(mb.group(3)) <= 1e-12 and mb.group(2) == "0", f"内部の面 {m.group(1)}、境界の半割面 {mb.group(3)} (対応の無い面 {mb.group(2)})")
else: und("2a.rSurfVect", "照合の行が無い")
mm = re.search(r"== 2b .*?\n(.*?)  rc=(\d+)", log, re.S); v("2b.float_vs_fp64", bool(mm) and mm.group(2) == "0", "rc=" + (mm.group(2) if mm else "なし"))
mm = re.search(r"== 2e .*?\n(.*?)  rc=(\d+)", log, re.S); v("2e.case39_3D_unchanged", bool(mm) and mm.group(2) == "0", "rc=" + (mm.group(2) if mm else "なし"))
def changed(a, b, allowed, extra_ok=()):
    fa, fb = h5py.File(a, "r"), h5py.File(b, "r"); nInt = int(fb["MESH"].attrs["nNormalPlanes"])
    names = []
    fa.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) else None)
    bad = []; ch = {}
    for n in names:
        if n not in fb: bad.append(f"{n} が新に無い"); continue
        x, y = np.asarray(fa[n][:]), np.asarray(fb[n][:])
        if x.shape != y.shape or x.dtype != y.dtype: bad.append(f"{n} の形・型"); continue
        diff = x.view(np.uint8) != y.view(np.uint8) if x.dtype.kind in "fiu" else (x != y)
        if np.any(diff):
            if n not in allowed: bad.append(n); continue
            k = 3 if x.size == 3 * fb["PLANES/surfArea"].shape[0] else 1
            idx = np.unique(np.flatnonzero(x != y) // k)
            ch[n] = (int(idx.size), int(idx.min()), float(np.max(np.abs(x.astype(np.float64) - y.astype(np.float64)))))
            if idx.min() < nInt: bad.append(f"{n} が内部の面 ({idx.min()} < {nInt}) で違う")
    newonly = []
    fb.visititems(lambda n, o: newonly.append(n) if isinstance(o, h5py.Dataset) and n not in fa else None)
    bad += [f"新だけにある {n}" for n in newonly if n not in extra_ok]
    return bad, ch, newonly
ALLOWED = {"PLANES/surfVect", "PLANES/surfArea", "PLANES/centCoords"}
for key, a, b, extra in (("2c.case45_changes", C45 / "_conv614/c45/new64/nozzle.h5", W / "c45/new64/nozzle.h5", ("PLANES/rSurfVect",)),
                         ("2d.case48_changes", C45 / "_conv614/c48/new64/m.h5", W / "c48/new64/m.h5", ())):
    try:
        bad, ch, newonly = changed(a, b, ALLOWED, extra)
        v(key, not bad, f"変わった (境界の半割面のみ許す): {ch}; 新だけ {newonly}; 許されない差 {bad}")
    except Exception as e: und(key, str(e))
# fg9.sh の run_0463〜0468 は化学種が移らず静止場が一様でなかった (fg9d.sh の冒頭) ので、取り直した run_0476〜0481 で判定する
print("== 3 デバイスの幾何の閉性 (FORGE_DIAG_HOOP_CLOSURE)")
def hoop(r):
    with h5py.File(C45 / r / "hoop.h5", "r") as h:
        nc = int(h.attrs["nCells"])   # A_planar・volume は nCells_all (周期・ゴースト込み) の長さ。判定は自前の CV (nCells) だけ
        d = {k: np.asarray(h["/cells/" + k][:nc], np.float64) for k in ("E_x", "E_y", "abs64_x", "abs64_y", "A_planar", "sum64_x", "sum64_y")}
        d["cc"] = np.asarray(h["/cells/cc64"][:], np.float64).reshape(-1, 3) if "/cells/cc64" in h else None
        d["attrs"] = {k: (h.attrs[k].item() if hasattr(h.attrs[k], "item") else h.attrs[k]) for k in h.attrs}
        d["ss"] = np.asarray(h["/faces/ss"][:], np.float64); d["nl"] = np.asarray(h["/faces/normlen_err"][:], np.float64)
    return d
def regions(n, cc):
    j = np.arange(n) % NJ; i = np.arange(n) // NJ; xr = cc[:n, 0] / RT
    return {"軸 j0": j == 0, "j 2〜8": (j >= 2) & (j <= 8), "収縮部の内部": (xr >= -5) & (xr < -1) & (j >= 20) & (j <= 60),
            "壁際 j117〜119": (j >= 117) & (j <= 119), "壁 j120": j == 120, "入口の角": (i == 0) & ((j <= 2) | (j >= 118))}
H = {}
try:
    for r in ("run_0476_hp_still64_k1", "run_0477_hp_still64_k0", "run_0478_hp_still32_k1", "run_0479_hp_still32_k1b", "run_0480_hp_still32_k0", "run_0481_hp_still32_k0b"): H[r] = hoop(r)
    k1 = H["run_0476_hp_still64_k1"]; n = k1["E_x"].size; R = regions(n, k1["cc"])
    ox = int(np.count_nonzero(np.abs(k1["sum64_x"]) > 100 * EPS * (k1["A_planar"] + k1["abs64_x"]))); oy = int(np.count_nonzero(np.abs(k1["sum64_y"] - k1["A_planar"]) > 100 * EPS * (k1["A_planar"] + k1["abs64_y"])))
    e28 = max(k1["E_x"][R["j 2〜8"]].max(), k1["E_y"][R["j 2〜8"]].max())
    fin = bool(np.all(np.isfinite(k1["ss"])) and np.all(k1["ss"] > 0))
    v("3.fp64_k1_closure", ox == 0 and oy == 0 and e28 <= 1e-10 and fin, f"超過 x {ox}・y {oy}、j 2〜8 の最大 {e28:.3e}、ss 有限・正 {fin}、法線の長さの誤差の最大 {np.max(np.abs(k1['nl'])):.2e}、ON {k1['attrs'].get('segment_rweight_active')} ({k1['attrs'].get('segment_rweight_reason')})")
    k0 = H["run_0477_hp_still64_k0"]
    print("   (記録) FP64 の領域ごとの最大 E_x / E_y  キー 1 | キー 0")
    for nm, m_ in R.items(): print(f"     {nm:10s} {k1['E_x'][m_].max():.2e} / {k1['E_y'][m_].max():.2e} | {k0['E_x'][m_].max():.2e} / {k0['E_y'][m_].max():.2e}")
    a1, a1b, a0, a0b = (H[r] for r in ("run_0478_hp_still32_k1", "run_0479_hp_still32_k1b", "run_0480_hp_still32_k0", "run_0481_hp_still32_k0b"))
    bad = []
    print("   float の領域ごと (最大 / RMS): キー 1 | キー 0 | 再実行の差")
    for nm, m_ in R.items():
        for c in ("E_x", "E_y"):
            mx1, mx0 = a1[c][m_].max(), a0[c][m_].max(); rm1, rm0 = np.sqrt(np.mean(a1[c][m_] ** 2)), np.sqrt(np.mean(a0[c][m_] ** 2))
            nz = max(abs(a1[c][m_].max() - a1b[c][m_].max()), abs(a0[c][m_].max() - a0b[c][m_].max()))
            nzr = max(abs(np.sqrt(np.mean(a1[c][m_] ** 2)) - np.sqrt(np.mean(a1b[c][m_] ** 2))), abs(rm0 - np.sqrt(np.mean(a0b[c][m_] ** 2))))
            ok = mx1 <= 1.1 * mx0 + 10 * nz and rm1 <= 1.1 * rm0 + 10 * nzr
            if not ok: bad.append(f"{nm} {c}")
            print(f"     {nm:10s} {c}: {mx1:.2e}/{rm1:.2e} | {mx0:.2e}/{rm0:.2e} | {nz:.1e}/{nzr:.1e} {'' if ok else '← 超過'}")
    v("3.float_nondegradation", not bad, f"超えた領域・成分 {bad}")
except Exception as e: und("3.closure", str(e))
print("== 4 一様な静止場 (自由な DOF = 軸 j0 と壁 j120 を除く)")
def precondition(r):   # 静止場の前提: res_0 の P・T・Y が一様で P = P0 (fg9.sh の誤りを再発させない)
    with h5py.File(C45 / r / "res_0.h5", "r") as h:
        V_ = h["VALUE"]; P = np.asarray(V_["P"][:], np.float64); T = np.asarray(V_["T"][:], np.float64)
        ys = [np.asarray(V_[k][:], np.float64) for k in V_ if re.fullmatch(r"Y\d+", k)]
        u = sum(np.abs(np.asarray(V_[k][:], np.float64)) for k in ("Ux", "Uy", "Uz"))
    rel = lambda a: (a.max() - a.min()) / abs(a.mean())
    msg = f"P {P.min():.9g}..{P.max():.9g}、T {T.min():.9g}..{T.max():.9g}、Y の幅の最大 {max((y.max() - y.min() for y in ys), default=np.nan):.3e}、|u| の最大 {u.max():.3e}"
    ok = bool(ys) and abs(P.mean() / P0 - 1) <= 1e-6 and rel(P) <= 1e-6 and rel(T) <= 1e-6 and max(y.max() - y.min() for y in ys) <= 1e-6 and u.max() == 0
    return ok, msg
def still(r):
    with h5py.File(C45 / r / "res_1.h5", "r") as h: q = np.asarray(h["VALUE/res_roUy"][:], np.float64)
    with h5py.File(C45 / r / "res_10.h5", "r") as h:
        V_ = h["VALUE"]; u = np.sqrt(sum(np.asarray(V_[k][:], np.float64) ** 2 for k in ("Ux", "Uy", "Uz")))
        bad = sum(int(np.count_nonzero(~np.isfinite(np.asarray(V_[k][:])))) for k in V_.keys())
    return q, u, bad
try:
    for r in H:
        ok, msg = precondition(r); print(f"   前提 {r}: {'OK' if ok else 'NG'} — {msg}")
        if not ok: raise ValueError(f"{r} の res_0 が一様な静止場でない ({msg})")
    S = {r: still(r) for r in H}
    Ap = H["run_0476_hp_still64_k1"]["A_planar"]; n = Ap.size; j = np.arange(n) % NJ; free = (j > 0) & (j < 120); R = regions(n, H["run_0476_hp_still64_k1"]["cc"])
    def stats(r, m_):
        q, u, bad = S[r]; z = np.abs(q / (P0 * Ap[:q.size]))[m_[:q.size] & free[:q.size]]
        return z.max(), np.sqrt(np.mean(z ** 2)), u[m_[:u.size] & free[:u.size]].max(), bad
    q28, _, _, b1 = stats("run_0476_hp_still64_k1", R["j 2〜8"]); _, _, umax1, _ = stats("run_0476_hp_still64_k1", np.ones(n, bool))
    v("4.fp64_k1_still", q28 <= 1e-9 and umax1 <= 2.6e-3 and b1 == 0, f"自由な DOF の j 2〜8 の |res_roUy|/(P·A) 最大 {q28:.3e} (≤ 1e-9)、10 step 後の最大 |u| {umax1:.3e} m/s (≤ 2.6e-3)、非有限 {b1}")
    q28k0, _, _, _ = stats("run_0477_hp_still64_k0", R["j 2〜8"]); _, _, umax0, _ = stats("run_0477_hp_still64_k0", np.ones(n, bool))
    print(f"   (記録) FP64 キー 0 (新しい格子): j 2〜8 の最大 {q28k0:.3e}、10 step 後の最大 |u| {umax0:.3e} m/s (旧の格子の run_0401 は 3.29e-3・0.257 m/s)")
    bad = []
    for nm, m_ in R.items():
        if nm == "軸 j0" or nm == "壁 j120": continue
        s1, s1b, s0, s0b = (stats(r, m_) for r in ("run_0478_hp_still32_k1", "run_0479_hp_still32_k1b", "run_0480_hp_still32_k0", "run_0481_hp_still32_k0b"))
        for idx, lab in ((0, "最大"), (1, "RMS"), (2, "|u|")):
            nz = max(abs(s1[idx] - s1b[idx]), abs(s0[idx] - s0b[idx]))
            ok = s1[idx] <= 1.1 * s0[idx] + 10 * nz
            if not ok: bad.append(f"{nm} {lab}")
            print(f"     float {nm:10s} {lab}: キー 1 {s1[idx]:.3e} | キー 0 {s0[idx]:.3e} | 再実行 {nz:.1e} {'' if ok else '← 超過'}")
    v("4.float_still_nondegradation", not bad and all(S[r][2] == 0 for r in S), f"超えた {bad}")
except Exception as e: und("4.still", str(e))
print("== 5 適用の分岐")
def faces(r):
    with h5py.File(C45 / r / "hoop.h5", "r") as h: return {k: np.asarray(h["/faces/" + k][:]) for k in ("sx", "sy", "sz", "ss")}, {k: h.attrs[k] for k in h.attrs}
try:
    for a, b, lab in (("run_0469_hp_br_floor_k1", "run_0470_hp_br_floor_k0", "axisRFloor 1e-4"), ("run_0471_hp_br_m1_k1", "run_0472_hp_br_m1_k0", "axisymMethod 1")):
        fa, aa = faces(a); fb, ab = faces(b)
        same = all(np.array_equal(fa[k].view(np.uint8), fb[k].view(np.uint8)) for k in fa)
        v(f"5.branch.{lab}", same, f"キー 1 と 0 で面の幾何がビット一致 {same}")
    def r1(r):
        with h5py.File(C45 / r / "res_1.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in ("res_ro", "res_roUx", "res_roUy", "res_roe", "ro", "roUx", "roUy", "roe")}
    A_, B_, Bb = r1("run_0474_hp_br_old"), r1("run_0473_hp_br_k0new"), r1("run_0475_hp_br_oldb")
    worst = 0.0
    for k in A_:
        s = np.max(np.abs(A_[k])); dn = np.max(np.abs(B_[k] - A_[k])) / s; dr = np.max(np.abs(Bb[k] - A_[k])) / s
        worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf))
    v("5.k0new_vs_old_same_h5", worst <= 10, f"新のキー 0 と旧のバイナリ (同じ新しい HDF5、1 step) の差 / 再実行の差の最大 {worst:.2f}")
except Exception as e: und("5.branch", str(e))
print("== 6 軸対称でないケースの回帰 (float、20 step)")
# case/48 の run_0051〜0053 (一様な初期値から) は新旧とも同じく step 9 で発散したので、収束場から始めた run_0054〜0056 (fg9c.sh) で判定する
for key, base, (o, nw, o2) in (("6.case48", C48, ("run_0055_hp_reg_old", "run_0054_hp_reg_new", "run_0056_hp_reg_oldb")), ("6.case56", C56, ("run_0094_hp_reg_old", "run_0093_hp_reg_new", "run_0095_hp_reg_oldb"))):
    try:
        def vals(r):
            with h5py.File(base / r / "res_20.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T")}
        A_, B_, A2 = vals(o), vals(nw), vals(o2); worst = 0.0
        for k in A_:
            s = np.max(np.abs(A_[k])); dn = np.max(np.abs(B_[k] - A_[k])) / s; dr = np.max(np.abs(A2[k] - A_[k])) / s
            if not (np.all(np.isfinite(B_[k]))): raise ValueError(f"{nw} の {k} に非有限")
            worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf))
        v(key, worst <= 10, f"新旧の差 / 再実行の差の最大 {worst:.2f}")
    except Exception as e: und(key, str(e))
json.dump(V, open(C45 / "fg9_verdict.json", "w"), ensure_ascii=False, indent=1, default=str)
bad = [k for k, x in V.items() if x["verdict"] != "PASS"]
print("== 総合:", "ALL PASS" if not bad else f"PASS でないもの {bad}")
