"""plan condensation-two-phase-default §5.1 #4g3 (G3-b 更新写像の収支) の集計と float 格納 / double 参照の A/B (最初の一手)。
  python3 g3b_judge.py <tp_update.h5> [<tp_update_2.h5> ...] [--length-scale 1000] [--csv out.csv]

入力は FORGE_DIAG_TP_UPDATE=<h5> で forge が書いた記録 (main.cpp runTpUpdateDiag; スロット配置は
solver_density_cuda/cuda_forge/twoPhaseUpdateDiag_d.cuh)。ON (run_0561 の設定) と OFF (run_0567 の設定) の h5 を並べて渡してよい。

恒等式 (Ω ごと・成分ごと、すべて double):
  Σ V (q_after − q_before) = Σ V δq_limited + Σ V C_round + Σ_a Σ V C_a + E_bookkeeping
  - q_before = 写し upd_start、q_after = 写し final (格納値)。
  - δq_limited と C_round は実効の commit (ON: TPC [総水分・液・Q]、OFF: SPC [総水分]・CMC [液・Q]) のもの。
    C_round = cand_f − cand_d (cand_d = before + δ を double で)、commit 内の床・vround は C_floor (C_a の 1 つ)。
  - C_a = 各補正操作がカーネル内で読んだ値と書いた値の差 (after − before)。ON の総水分は、化学種の commit (SPC; δ・C_round・床を含む全差)
    と hold (水の戻し; 写し hold − sp_commit) を明示し、両者の和が 0 になるはずの組として別行にする (二重計上しない)。
  - E_bookkeeping = 左辺 − 右辺。操作をすべて記録していれば double の丸めだけになる。0 から離れていれば記録点の欠落
    (操作と操作の間で格納値を書き換えた未記録の処理) を示す。境目の照合 (gap) で場所を出す。
  - ρv = ρY_w − ρg は各項を double で差し引いて作る。
A/B (codex 2026-10-04 の判別): A = 本番の float 格納のまま C_round を計上しない未閉鎖 (左辺 − Σδ − ΣC_a)、
  B = commit の候補を double で格納したときの未閉鎖 (A − ΣC_round)。A − B = ΣC_round。B が double の丸め程度なら
  「A の未閉鎖は格納丸めで説明される」、B が残れば補正項・記録点の欠落を追う。
  C_round の大きさは IEEE の最近接丸めの上界 (u = 2⁻²⁴、非正規化の半刻み 2⁻¹⁵⁰) と節点ごとに照合し、超えた節点数を出す (参考)。
合否の閾値は置かない (plan: この A/B の後に決める)。|E_bookkeeping| / Σ|項| と A/B の数値を表で出すだけ。

Ω: all、x 帯 [入口, 17)・[17, 35)・[35, 70)・[70, 出口] mm、x 35–70 mm の壁距離 < 0.4 mm と ≥ 0.4 mm、参考に入口ピン節点。
長さは h5 の ccx・wall_dist に --length-scale (既定 1000 = m → mm) を掛けた値で切る (x の範囲を最初に表示するので単位を確かめること)。
"""
import argparse, csv, sys
import h5py
import numpy as np

U = 2.0 ** -24          # float32 の単位丸め (最近接)
SUBHALF = 2.0 ** -150   # 非正規化数の半刻み
COMPS = ["w", "g", "Q2", "Q1", "Q0"]


class Rec:
    def __init__(self, path, scale):
        f = h5py.File(path, "r")
        self.path = path
        self.names = [s.decode() if isinstance(s, bytes) else s for s in f["/slot_names"][...]]
        self.idx = {nm: i for i, nm in enumerate(self.names)}
        self.upd = f["/update/slots"][...]
        self.pre = f["/pre/slots"][...]
        self.labels = [s.decode() if isinstance(s, bytes) else s for s in f["/snap/labels"][...]]
        self.snap = f["/snap/data"][...].astype(np.float64)   # [ns, 5, n]
        self.V = f["/node/volume"][...].astype(np.float64)
        self.x = f["/node/ccx"][...].astype(np.float64) * scale
        self.wd = f["/node/wall_dist"][...].astype(np.float64) * scale if "/node/wall_dist" in f else None
        self.pin = (f["/node/scalarDirichletPin"][...] == 1.0) if "/node/scalarDirichletPin" in f else np.zeros_like(self.V, bool)
        self.on = int(f.attrs["twophase_active"]) == 1
        self.attrs = {k: f.attrs[k] for k in f.attrs.keys()}
        self.n = self.V.size
        f.close()

    def s(self, phase, name):
        buf = self.upd if phase == "upd" else self.pre
        return buf[self.idx[name]]

    def snapshot(self, label, c):
        if label not in self.labels:
            return None
        return self.snap[self.labels.index(label), c]


def opn(rec, short):
    for nm in rec.attrs["ops"]:
        nm = nm.decode() if isinstance(nm, bytes) else nm
        if nm.startswith(short + "_"):
            return nm
    raise KeyError(short)


def ba(rec, phase, short, c):
    """操作 short の成分 c の (before, after)。経路に無ければ None。"""
    nm = opn(rec, short)
    b = rec.s(phase, f"{nm}/{COMPS[c]}/before")
    a = rec.s(phase, f"{nm}/{COMPS[c]}/after")
    if not np.isfinite(b).any():
        return None
    return b, a


def commit_x(rec, short, c, x):
    return rec.s("upd", f"{opn(rec, short)}/{COMPS[c]}/{x}")


def z(a):
    return np.where(np.isfinite(a), a, 0.0)


# 経路 (操作の順)。HOLD は写しから作る。
def path_ops(rec, c):
    if c == 0:
        return ["SPC", "HOLD", "TPC", "RNF", "RNS"] if rec.on else ["SPC", "RNF", "RNS"]
    tail = ["PFL", "RZL", "RZU", "RZP", "RZR"]
    return (["TPC", "RNS"] + tail) if rec.on else (["CMC", "ARH", "LIM"] + tail)


def node_terms(rec, c):
    """成分 c の節点ごとの項 (dict 名前 → 配列) と、境目の照合結果、C_round の上界。"""
    q0 = rec.snapshot("upd_start", c); q1 = rec.snapshot("final", c)
    T = {"LHS": q1 - q0}
    eff = "TPC" if rec.on else ("SPC" if c == 0 else "CMC")
    gaps = []; prev_after = q0; prev_name = "upd_start"
    bound = None; cr = None
    for op in path_ops(rec, c):
        if op == "HOLD":
            b = rec.snapshot("sp_commit", c); a = rec.snapshot("hold", c)
            if b is None or a is None:
                continue
            T["C_hold (undo species commit)"] = a - b
        else:
            r = ba(rec, "upd", op, c)
            if r is None:
                continue
            b, a = r
            pres = np.isfinite(b) & np.isfinite(a)
            if op in ("SPC", "TPC", "CMC"):
                dd = z(commit_x(rec, op, c, "delta_d")); cf = commit_x(rec, op, c, "cand_f"); cd = commit_x(rec, op, c, "cand_d")
                crr = z(cf - cd); flo = z(a - cf)
                lost = z(commit_x(rec, op, c, "lost"))
                if op == eff:
                    T["delta_limited"] = dd; T["C_round"] = crr; T["C_floor (commit)"] = flo; T["lost_to_storage"] = lost
                    cr = crr
                    # IEEE 上界: 格納 1 回 + (TPC は積・内側の和の丸め; 本番が融合していれば過大側)
                    bnd = U * np.abs(z(cd)) + SUBHALF
                    if op == "TPC":
                        if c == 0:
                            dg = z(commit_x(rec, "TPC", 1, "delta_d")); dv = dd - dg
                            bnd = bnd + U * np.abs(dd) + U * (np.abs(dv) + np.abs(dg)) + 2 * SUBHALF
                        else:
                            bnd = bnd + U * np.abs(dd) + SUBHALF
                    bound = bnd
                else:   # ON の水の化学種 commit (hold で戻る): 全差を 1 行にまとめる
                    T["C_species_commit (undone by hold)"] = z(a - b)
            else:
                key = {"RNF": "C_renorm_negfloor", "RNS": "C_renorm_scale", "PFL": "C_passive_floor", "ARH": "C_add_rho_term",
                       "LIM": "C_limit_increment", "RZL": "C_realiz_lower", "RZU": "C_realiz_upper", "RZP": "C_projection",
                       "RZR": "C_droplet_removal (physical)"}[op]
                T[key] = z(np.where(pres, a - b, 0.0))
        # 境目: この操作が読んだ値 = 直前の操作が書いた値 (ビット一致) か
        if op == "HOLD":
            b_here = rec.snapshot("sp_commit", c); a_here = rec.snapshot("hold", c)
        else:
            b_here, a_here = r
        m = np.isfinite(b_here)
        nmis = int(np.count_nonzero(m & (b_here != prev_after)))
        gaps.append((f"{prev_name} -> {op}", nmis, float(np.sum(rec.V[m] * np.abs(b_here[m] - prev_after[m])))))
        prev_after = np.where(np.isfinite(a_here), a_here, prev_after); prev_name = op
    m = np.isfinite(q1)
    gaps.append((f"{prev_name} -> final", int(np.count_nonzero(q1 != prev_after)), float(np.sum(rec.V * np.abs(q1 - prev_after)))))
    T["B_boundary_overwrite"] = np.zeros_like(q0)   # 更新の中に境界の上書きは無い (入口ピンは前処理; 下の前処理収支)
    return T, gaps, bound, cr


def regions(rec):
    x = rec.x
    R = [("all", np.ones(rec.n, bool)),
         ("x[in,17)", x < 17.0), ("x[17,35)", (x >= 17.0) & (x < 35.0)),
         ("x[35,70)", (x >= 35.0) & (x < 70.0)), ("x[70,ex]", x >= 70.0)]
    if rec.wd is not None:
        band = (x >= 35.0) & (x < 70.0)
        R += [("x35-70 wd<0.4", band & (rec.wd < 0.4)), ("x35-70 wd>=0.4", band & (rec.wd >= 0.4))]
    R += [("pinned (ref)", rec.pin)]
    return R


def fmt(v):
    return f"{v:+.4e}"


def update_budget(rec, rows_csv):
    tag = "ON" if rec.on else "OFF"
    terms = {}; gaps = {}; bounds = {}; crs = {}
    for c in range(5):
        terms[c], gaps[c], bounds[c], crs[c] = node_terms(rec, c)
    # ρv = w − g (項名ごとに double で差し引く)
    names_all = []
    for c in range(5):
        for k in terms[c]:
            if k not in names_all:
                names_all.append(k)
    tv = {k: terms[0].get(k, 0.0) - terms[1].get(k, 0.0) for k in names_all}
    cols = COMPS + ["v"]
    getT = lambda ci, k: (tv.get(k, 0.0) if ci == 5 else terms[ci].get(k, None))
    print(f"\n=== {rec.path}  [{tag}]  update-map budget (one outer implicit update) ===")
    for c in range(5):
        bad = [(g, n_, s_) for (g, n_, s_) in gaps[c] if n_ > 0]
        print(f"  boundary check {COMPS[c]:2s}: " + ("all op boundaries bit-identical" if not bad else
              "; ".join(f"{g}: {n_} nodes, sum V|d| {s_:.3e}" for g, n_, s_ in bad)))
    order = ["delta_limited", "C_round", "C_floor (commit)", "C_species_commit (undone by hold)", "C_hold (undo species commit)",
             "C_renorm_negfloor", "C_renorm_scale", "C_passive_floor", "C_add_rho_term", "C_limit_increment",
             "C_realiz_lower", "C_realiz_upper", "C_projection", "C_droplet_removal (physical)", "B_boundary_overwrite"]
    for rname, mask in regions(rec):
        if not mask.any():
            print(f"\n  -- Omega = {rname}: no nodes"); continue
        V = rec.V[mask]
        print(f"\n  -- Omega = {rname} ({int(mask.sum())} nodes); signed sum V*term [abs sum V*|term|] --")
        print("  " + f"{'term':38s}" + "".join(f"{c:>26s}" for c in cols))
        def S(ci, k, absol=False):
            t = getT(ci, k)
            if t is None or np.isscalar(t):
                return 0.0
            t = t[mask]
            return float(np.sum(V * (np.abs(t) if absol else t)))
        lhs = [S(ci, "LHS") for ci in range(6)]
        print("  " + f"{'LHS sum V(q_after-q_before)':38s}" + "".join(f"{fmt(v):>26s}" for v in lhs))
        rhs = [0.0]*6; absterms = [0.0]*6; dsum = [0.0]*6; crsum = [0.0]*6; casum = [0.0]*6
        for k in order:
            vals = [S(ci, k) for ci in range(6)]; avals = [S(ci, k, True) for ci in range(6)]
            if all(v == 0.0 for v in vals) and all(v == 0.0 for v in avals) and k != "B_boundary_overwrite":
                continue
            print("  " + f"{k:38s}" + "".join(f"{fmt(v)+' ['+f'{a:.2e}'+']':>26s}" for v, a in zip(vals, avals)))
            for ci in range(6):
                rhs[ci] += vals[ci]; absterms[ci] += abs(vals[ci])
                if k == "delta_limited": dsum[ci] += vals[ci]
                elif k == "C_round": crsum[ci] += vals[ci]
                else: casum[ci] += vals[ci]
            for ci in range(6):
                rows_csv.append([rec.path, tag, rname, cols[ci], k, vals[ci], avals[ci]])
        E = [lhs[ci] - rhs[ci] for ci in range(6)]
        lost = [S(ci, "lost_to_storage") for ci in range(6)]
        print("  " + f"{'E_bookkeeping = LHS - sum(terms)':38s}" + "".join(f"{fmt(v):>26s}" for v in E))
        print("  " + f"{'|E| / sum|terms|':38s}" + "".join(f"{(abs(e)/a if a > 0 else float('nan')):>26.3e}" for e, a in zip(E, absterms)))
        print("  " + f"{'increment lost to storage (sum V)':38s}" + "".join(f"{fmt(v):>26s}" for v in lost))
        nlost = [int(np.count_nonzero(terms[ci]["lost_to_storage"][mask] != 0)) if ci < 5 and "lost_to_storage" in terms[ci] else 0 for ci in range(5)]
        print("  " + f"{'  nodes with lost increment':38s}" + "".join(f"{v:>26d}" for v in nlost))
        # A/B
        A = [lhs[ci] - dsum[ci] - casum[ci] for ci in range(6)]
        B = [A[ci] - crsum[ci] for ci in range(6)]
        print("  A/B (A = float storage, C_round not booked; B = candidates stored in double):")
        print("  " + f"{'  A unclosed = LHS - sum delta - sum C_a':38s}" + "".join(f"{fmt(v):>26s}" for v in A))
        print("  " + f"{'  B unclosed = A - sum C_round':38s}" + "".join(f"{fmt(v):>26s}" for v in B))
        print("  " + f"{'  A - B (= sum V C_round)':38s}" + "".join(f"{fmt(a - b):>26s}" for a, b in zip(A, B)))
        print("  " + f"{'  |B| / |A|':38s}" + "".join(f"{(abs(b)/abs(a) if a != 0 else float('nan')):>26.3e}" for a, b in zip(A, B)))
        # IEEE 上界との照合 (参考)
        vb = []; vc = []
        for ci in range(5):
            if bounds[ci] is None:
                vb.append("n/a"); vc.append("n/a"); continue
            bm = bounds[ci][mask]; cm = np.abs(crs[ci][mask])
            vb.append(f"{int(np.count_nonzero(cm > bm))}"); vc.append(f"{float(np.sum(V*bm)):.3e}")
        print("  " + f"{'  C_round nodes above IEEE bound':38s}" + "".join(f"{v:>26s}" for v in vb))
        print("  " + f"{'  sum V*IEEE bound (ref)':38s}" + "".join(f"{v:>26s}" for v in vc))
        for ci in range(6):
            rows_csv.append([rec.path, tag, rname, cols[ci], "LHS", lhs[ci], ""])
            rows_csv.append([rec.path, tag, rname, cols[ci], "E_bookkeeping", E[ci], ""])
            rows_csv.append([rec.path, tag, rname, cols[ci], "A_unclosed", A[ci], ""])
            rows_csv.append([rec.path, tag, rname, cols[ci], "B_unclosed", B[ci], ""])


def pre_budget(rec, rows_csv):
    tag = "ON" if rec.on else "OFF"
    print(f"\n=== {rec.path}  [{tag}]  pre-part budget (assembly before the update) ===")
    print("  snapshots: " + " ".join(rec.labels))
    for rname, mask in regions(rec):
        if not mask.any():
            continue
        V = rec.V[mask]
        print(f"  -- Omega = {rname} --")
        print("  " + f"{'term':42s}" + "".join(f"{c:>14s}" for c in COMPS))
        out = {}
        for c in range(5):
            i0 = rec.snapshot("init", c); ps = rec.snapshot("pre_state", c); pe = rec.snapshot("pre_end", c); po = rec.snapshot("post_end", c)
            d = {"LHS sum V(pre_end - init)": pe - i0}
            rz = np.zeros_like(i0)
            for op, key in [("RZL", "realiz_lower"), ("RZU", "realiz_upper"), ("RZP", "projection"), ("RZR", "droplet_removal (phys)")]:
                r = ba(rec, "pre", op, c)
                t = np.zeros_like(i0) if r is None else z(r[1] - r[0])
                d[key] = t; rz += t
            d["other state ops (pre_state-init-RZ)"] = (ps - i0) - rz
            d["boundary overwrite (pre_end-pre_state)"] = pe - ps
            d["post-part change (post_end-pre_end)"] = po - pe
            out[c] = d
        for k in out[0]:
            vals = [float(np.sum(V * out[c][k][mask])) for c in range(5)]
            print("  " + f"{k:42s}" + "".join(f"{fmt(v):>14s}" for v in vals))
            for c in range(5):
                rows_csv.append([rec.path, tag, rname, COMPS[c], "pre:" + k, vals[c], ""])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5", nargs="+")
    ap.add_argument("--length-scale", type=float, default=1000.0, help="ccx・wall_dist に掛ける倍率 (既定 1000 = m → mm)")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    rows = []
    for p in a.h5:
        rec = Rec(p, a.length_scale)
        print(f"# {p}: {rec.n} nodes, twophase_active={int(rec.on)}, x range [{rec.x.min():.3f}, {rec.x.max():.3f}] (scaled; mm expected)"
              + ("" if rec.wd is None else f", wall_dist max {rec.wd.max():.3f}"))
        for k in ("cfl_pseudo", "implicitRelax", "nStepInner", "condTwoPhaseSolver", "condTwoPhaseNonnegLimit", "condTwoPhaseRelax", "input_value_file"):
            if k in rec.attrs:
                print(f"#   {k} = {rec.attrs[k]}")
        for l in rec.attrs.get("summary", []):
            print("#   " + (l.decode() if isinstance(l, bytes) else str(l)))
        pre_budget(rec, rows)
        update_budget(rec, rows)
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["file", "arm", "omega", "component", "term", "signed_sumV", "abs_sumV"]); w.writerows(rows)
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    sys.exit(main())
