"""plan architecture-float-state-double-geometry §6.14 (事前登録) の集計。fg8.sh の後に AWS の case/45 で回す。結果は fg8_verdict.json。
欠損・非有限は判定不能 (= PASS でない)。"""
import json, re, h5py, numpy as np
from pathlib import Path
C45 = Path("/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"); C48 = C45.parent / "48.flat_plate_cooled_m4"; C56 = C45.parent / "56.gap_tp1187"
V = {}
def v(key, ok, det): V[key] = {"verdict": "PASS" if ok else "FAIL", "detail": det}; print(f"  [{'PASS' if ok else 'FAIL'}] {key}: {det}")
log = (C45 / "fg8.log").read_text().split("== 開始")[-1]
print("== 1・2 変換の比較 (compare_mesh_h5.py の終了コード)")
for cs in ("c45", "c48", "c39"):
    conv = re.findall(rf"conv {cs} (\w+) rc=(\d+)", log)
    for n, tag in ((1, "新 FP64 対 旧 FP64 (全件)"), (2, "新 float 対 新 FP64 (--geometry-only --require-float64)")):
        m = re.search(rf"== {n} {cs} .*?\n(.*?)  rc=(\d+)", log, re.S)
        if m is None: v(f"{n}.{cs}", False, "結果の行が無い (判定不能)"); continue
        summ = [l for l in m.group(1).splitlines() if "summary" in l or "VERDICT" in l]
        v(f"{n}.{cs}", m.group(2) == "0" and all(rc == "0" for _, rc in conv), f"{tag}: rc={m.group(2)}、変換 {conv}; " + " / ".join(s.strip() for s in summ))
print("== 3 警告")
warn = {k: bool(re.search(rf"== 3 警告 {re.escape(k)}: .*float32", log)) for k in ("run_0048_fg6_dump32", "run_0049_fg6_dump64", "run_0050_fg6_warn_mixed", "case/56 (float32 の格子)")}
exp = {"run_0048_fg6_dump32": False, "run_0049_fg6_dump64": False, "run_0050_fg6_warn_mixed": True, "case/56 (float32 の格子)": True}
v("3.warning", warn == exp, f"出た: {warn} (期待 {exp})")
print("== 4 新 float のソルバが新しい double の格子を読む (段 ② のダンプ)")
try:
    f32 = h5py.File(C48 / "run_0048_fg6_dump32/stage2.h5", "r"); f64 = h5py.File(C48 / "run_0049_fg6_dump64/stage2.h5", "r")
    g = (int(f32.attrs["has_geom64"]), int(f64.attrs["has_geom64"]), int(f32.attrs["flow_float_bytes"]), int(f64.attrs["flow_float_bytes"]))
    c32 = np.asarray(f32["/lsq/new/cInt"][:], np.float64).reshape(-1, 3); c64 = np.asarray(f64["/lsq/new/cInt"][:], np.float64).reshape(-1, 3)
    if not (np.all(np.isfinite(c32)) and np.all(np.isfinite(c64))): raise ValueError("cInt に非有限")
    n64 = np.linalg.norm(c64, axis=1); nz = n64 > 0
    rv = float(np.max(np.linalg.norm(c32 - c64, axis=1)[nz] / n64[nz]))
    irep = sum(int(np.count_nonzero(np.asarray(f32[f"/wall/{i}/irep"][:]) != np.asarray(f64[f"/wall/{i}/irep"][:]))) for i in f64["/wall"].keys())
    d1 = max(float(np.nanmax(np.abs(np.asarray(f32[f"/d1d2/{i}/d1"][:]) - np.asarray(f64[f"/d1d2/{i}/d1"][:])) / np.abs(np.asarray(f64[f"/d1d2/{i}/d1"][:])))) for i in f64["/d1d2"].keys())
    v("4.double_mesh_read", g[0] == 1 and g[1] == 1 and g[2] == 4 and g[3] == 8 and rv <= 1e-6 and irep == 0 and d1 <= 1e-6,
      f"has_geom64 float/FP64 {g[0]}/{g[1]} (flow_float {g[2]}/{g[3]} B)、cInt のベクトルの相対差の最大 {rv:.3g}、irep の不一致 {irep}、d1 の相対差の最大 {d1:.3g}")
except Exception as e:
    v("4.double_mesh_read", False, f"判定不能 ({e})")
print("== 5 case/56 の 20 step の互換性")
try:
    def vals(r):
        with h5py.File(C56 / r / "res_20.h5", "r") as h:
            out = {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T")}
        for k, a in out.items():
            if not np.all(np.isfinite(a)): raise ValueError(f"{r} の {k} に非有限")
        return out
    A, B, A2 = vals("run_0091_fg6_id32_old"), vals("run_0090_fg6_id32_new"), vals("run_0092_fg6_id32_old2")
    worst = 0.0; parts = []
    for k in A:
        s = np.max(np.abs(A[k])); dn = np.max(np.abs(B[k] - A[k])) / s; dr = np.max(np.abs(A2[k] - A[k])) / s
        worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf)); parts.append(f"{k} {dn:.2e}/{dr:.2e}")
    v("5.compat.case56", worst <= 10, f"新旧の差 / 再実行の差の最大比 {worst:.2f}; " + "; ".join(parts))
except Exception as e:
    v("5.compat.case56", False, f"判定不能 ({e})")
json.dump(V, open(C45 / "fg8_verdict.json", "w"), ensure_ascii=False, indent=1)
bad = [k for k, x in V.items() if x["verdict"] != "PASS"]
print("== 総合:", "ALL PASS" if not bad else f"PASS でないもの {bad}")
