"""plan architecture-float-state-double-geometry §6.62 (試験 P、事前登録・2026-10-13 の資格の訂正の後) の、1 出力の取り出し。
res_<N>.h5 から局所の残差と物理性を 1 行の JSON にして <run>/tp_series.jsonl に足す。tp_arm.sh が出力ごとに呼び、取り出せたら場を消す。

局所の指標 (段 1 の打ち切りの規則): 区間 x/r_t ∈ [−8, −2) の全内部節点の Σ|res_q| (q = ro・roUx・roUy・roe・roK・roOmega)。
  - 取り出し方は forge-25 の cg_region_res.py (origin/feature/nozzle-core-grid 85bcd83e) を参考にしたが、量は登録どおり Σ|res_q|
    (cg_region_res.py の体積の重み・RMS はかけない)。
  - 内部節点: 節点を構造の並び (ni, nj) = (4467, 121) (G0c、cold_xcheck.mesh_info の EXTRA_GRIDS) に並べ、外周を境界の節点として除く:
      i = 0 (入口)・i = ni − 1 (出口)・j = 0 (軸)・j = nj − 1 (壁)。
    向き (j = 0 が軸、j = nj − 1 が壁、i とともに x が増える) は cold_xcheck.reduce_fields と同じで、毎回座標で確かめる (違えば止める)。
    区間 [−8, −2) に入口・出口の列は入らないが、同じ規則で除く。
  - x/r_t は節点ごとの x を cold_xcheck.mesh_info の S (r_t) で割った値 (壁法線の線は一定 x ではないので、列でなく節点で区間に入れる)。
  - 区間と除き方はコードの定数 (X_LO・X_HI・外周)。起動の前に `--fix <入力 h5>` で節点の集合 (数・添字の sha256) を
    _band_ab/cold_pair/tp_region_set.json に固定し、各出力ではその run の nozzle.h5 から作り直した集合が固定したものと同じことを確かめる
    (違えば取り出さない)。
物理性 (資格 Q6): 全節点の ρ・P・T の非有限の数・最小・0 以下の数 (cg_phys.py の量に 0 以下の数を足した)、
  Y0 の最大 − 最小 (VALUE/Y0、無ければ roY0/ro)、VALUE の全浮動小数点の量の非有限の数 (記録)。
非有限の残差は和を null にして数を別に書く (判定器が比の判定から分ける)。
usage:
  python3 tp_region.py --fix <入力 h5>     (起動の前に 1 回。既に固定した集合と違えば失敗)
  python3 tp_region.py <run> <N>           (各出力。同じ N の行が既にあれば何もしない。取り出せなければ 0 以外で終わり、行を書かない)
"""
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402

OUTD = HERE / "_band_ab" / "cold_pair"
SETF = OUTD / "tp_region_set.json"
X_LO, X_HI = -8.0, -2.0                 # 区間 [X_LO, X_HI) (r_t 単位、§6.62)
GRID = (4467, 121)                      # G0c (ni, nj)
RES_Q = ("ro", "roUx", "roUy", "roe", "roK", "roOmega")
AXIS_TOL = 1e-6                         # 軸の行の |r|/r_t の上限 (向きの確認)


def region_mask(xy):
    """節点座標 (n, ≥2) → 区間の内部節点の印 (n,)、S。向きが想定と違えば止める。"""
    ni, nj, S = XC.mesh_info(len(xy))
    if (ni, nj) != GRID:
        raise SystemExit(f"格子が G0c (ni {GRID[0]}・nj {GRID[1]}) でない: ({ni}, {nj})")
    x = xy[:, 0].reshape(ni, nj) / S
    r = xy[:, 1].reshape(ni, nj) / S
    if not np.max(np.abs(r[:, 0])) <= AXIS_TOL:
        raise SystemExit(f"j = 0 が軸 (r = 0) でない (max |r|/r_t {np.max(np.abs(r[:, 0])):.3e})")
    if not np.all(r[:, -1] > r[:, -2]):
        raise SystemExit("j = nj − 1 が壁 (各列で r が最大) でない")
    if not np.all(np.diff(x[:, -1]) > 0):
        raise SystemExit("壁の x が i とともに増えない")
    inner = np.zeros((ni, nj), bool)
    inner[1:-1, 1:-1] = True            # 外周 (入口 i=0・出口 i=ni−1・軸 j=0・壁 j=nj−1) を除く
    m = inner & (x >= X_LO) & (x < X_HI)
    return m.ravel(), x


def set_id(m, x):
    idx = np.flatnonzero(m).astype("<i8")
    xs = x.ravel()[idx]
    ii = idx // GRID[1]; jj = idx % GRID[1]
    return {"n": int(idx.size), "sha16": hashlib.sha256(idx.tobytes()).hexdigest()[:16], "x_lo": X_LO, "x_hi": X_HI,
            "grid": list(GRID), "x_range": [float(xs.min()), float(xs.max())] if idx.size else None,
            "i_range": [int(ii.min()), int(ii.max())] if idx.size else None, "j_range": [int(jj.min()), int(jj.max())] if idx.size else None,
            "boundary_rule": "外周 (i=0 入口・i=ni−1 出口・j=0 軸・j=nj−1 壁) を除く"}


def coords(p):
    with h5py.File(p, "r") as h:
        return np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)


def fix(inp):
    m, x = region_mask(coords(inp))
    sid = set_id(m, x)
    if sid["n"] == 0:
        raise SystemExit("区間の内部節点が 0")
    sid["source"] = str(inp)
    if SETF.exists():
        old = json.loads(SETF.read_text())
        if (old.get("n"), old.get("sha16")) != (sid["n"], sid["sha16"]):
            raise SystemExit(f"既に固定した集合 ({old.get('n')}, {old.get('sha16')}) と違う ({sid['n']}, {sid['sha16']})")
        print(f"[tp_region] 同じ集合が既に固定されている: {sid['n']} 節点 sha {sid['sha16']}")
        return
    OUTD.mkdir(parents=True, exist_ok=True)
    SETF.write_text(json.dumps(sid, ensure_ascii=False, indent=1))
    print(f"[tp_region] 集合を固定した: {sid['n']} 節点 sha {sid['sha16']} x/r_t {sid['x_range']} i {sid['i_range']} j {sid['j_range']} → {SETF}")


def finite_min(a):
    f = a[np.isfinite(a)]
    return float(f.min()) if f.size else None


def extract(run, n):
    d = HERE / run
    out = d / "tp_series.jsonl"
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip() and json.loads(line).get("step") == n:
                print(f"[tp_region] {run} step {n}: 既に行がある (何もしない)")
                return
    if not SETF.exists():
        raise SystemExit(f"集合が固定されていない ({SETF} が無い。起動の前に --fix を回す)")
    fixed = json.loads(SETF.read_text())
    m, x = region_mask(coords(d / "nozzle.h5"))
    sid = set_id(m, x)
    if (sid["n"], sid["sha16"]) != (fixed["n"], fixed["sha16"]):
        raise SystemExit(f"節点の集合が固定したものと違う ({sid['n']}, {sid['sha16']} 対 {fixed['n']}, {fixed['sha16']})")
    rec = {"run": run, "step": n, "set_n": sid["n"], "set_sha16": sid["sha16"]}
    with h5py.File(d / f"res_{n}.h5", "r") as h:
        V = h["VALUE"]
        need = ["ro", "P", "T"] + ([] if "Y0" in V else ["roY0"])
        miss = [k for k in need if k not in V]
        if miss:
            raise SystemExit(f"必要な量が無い: {miss}")
        if len(V["ro"]) != len(m):
            raise SystemExit(f"節点の数が格子と違う ({len(V['ro'])} 対 {len(m)})")
        mres = [q for q in RES_Q if "res_" + q not in V]
        if mres and n > 0:
            raise SystemExit(f"残差の量が無い: {['res_' + q for q in mres]}")
        rec["sum_abs_res"] = {}; rec["nonfinite_res"] = {}
        for q in RES_Q:
            if q in mres:                 # step 0 だけ許す (判定の窓は step 15,000 以上)
                rec["sum_abs_res"][q] = None; rec["nonfinite_res"][q] = None
                continue
            a = np.asarray(V["res_" + q][:], float)[m]
            nf = int(np.count_nonzero(~np.isfinite(a)))
            rec["sum_abs_res"][q] = float(np.sum(np.abs(a))) if nf == 0 else None
            rec["nonfinite_res"][q] = nf
        if mres:
            rec["missing_res"] = mres
        for k in ("ro", "P", "T"):
            a = np.asarray(V[k][:], float)
            fin = np.isfinite(a)
            rec[f"nonfinite_{k}"] = int(np.count_nonzero(~fin))
            rec[f"min_{k}"] = finite_min(a)
            rec[f"nonpos_{k}"] = int(np.count_nonzero(a[fin] <= 0))
        if "Y0" in V:
            y = np.asarray(V["Y0"][:], float); rec["Y0_src"] = "VALUE/Y0"
        else:
            y = np.asarray(V["roY0"][:], float) / np.asarray(V["ro"][:], float); rec["Y0_src"] = "roY0/ro"
        rec["Y0_nonfinite"] = int(np.count_nonzero(~np.isfinite(y)))
        rec["Y0_ptp"] = float(np.ptp(y)) if rec["Y0_nonfinite"] == 0 else None
        rec["Y0_min_max"] = [finite_min(y), (float(np.max(y[np.isfinite(y)])) if np.isfinite(y).any() else None)]
        rec["nonfinite_all_values"] = int(sum(np.count_nonzero(~np.isfinite(np.asarray(V[k][:], float)))
                                              for k in V if V[k].dtype.kind == "f"))
    with open(out, "a") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    s = rec["sum_abs_res"]
    print(f"[tp_region] {run} step {n}: Σ|res_ro| {s['ro']} Σ|res_roe| {s['roe']} 非有限 ρPT "
          f"{rec['nonfinite_ro'] + rec['nonfinite_P'] + rec['nonfinite_T']} Y0 ptp {rec['Y0_ptp']}")


def main(argv):
    if len(argv) == 2 and argv[0] == "--fix":
        fix(Path(argv[1]) if Path(argv[1]).is_absolute() else HERE / argv[1])
    elif len(argv) == 2:
        extract(argv[0], int(argv[1]))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
