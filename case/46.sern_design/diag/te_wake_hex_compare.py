#!/usr/bin/env python3
"""後縁下流の局所変形の 0 step 比較: A (te_wake_blend_H 0) と B (1.0) の最終 node 格子で、対応するヘキサごとの
最小スケール済みヤコビアンと AR を比べる (plan tooling-sern-te-wake-grid §5.1 #4、codex diagnose 2026-10-11 g4-admission-jacobian)。
読むだけ。変形領域 = 動いた節点を 1 個以上含むヘキサ。
使い方: te_wake_hex_compare.py <A の run> <B の run> [--thr 0.65] [--ar-max 5000] [--json out.json]

変位の局所化 (plan §5.1 #4c、事前登録 2026-10-11、codex diagnose g4-grid-redesign) の判定には **A = 0 格子 (参照、run_1080 相当)・
B = 局所化 (`mesh3d.te_wake_mode: local`)** で使う (判定 (1) は 0 格子との対応比較)。現方式 (band) の対照は A = 0 格子・B = 現方式で
別に回す (run_1081 の記録 `notes/investigations/2026-10-11-g4-admission/HEX_COMPARE.json`)。この道具が出す #4c の行:
  (1) 新しく 0.65 未満 0・既存の 0.65 未満の悪化 0・新しい AR 超過 (A ≤ --ar-max・B > --ar-max) 0 (全域の件数も残す)
  (2) B の変形領域の最小スケール済みヤコビアン ≥ --thr・全域で非正/非有限の頂点ヤコビアン 0 (向きは A の多数派)
  (4 の一部) 支持領域の外の座標が A とビット一致。支持領域 = B の info の変形区間の station (te_wake_i_first..te_wake_i_last) の
      主ブロックで、旧中間線 (A の j = jm) からの鉛直距離 < TE_WAKE_LOCAL_D_ZERO_H (0.30 H、float32 の丸めの余裕 1e-6 H)。
      B の info の te_wake_mode_version が現行のメッシャの定数と違えば判定不能
  最後に `TE_WAKE_LOCAL VERDICT: PASS|FAIL|UNDECIDABLE` (終了コードは従来どおり 0、接続・節点数の不一致だけ 2)。
(3) の後縁・復帰区間の折れ・層厚・固定座標と境界・skew と (4) の check_mesh_quality・check_dual_closure は
te_wake_grid_check.py --admission の行で見る (同じ A・B で回す)。
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import te_wake_grid_check as T  # noqa: E402

EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
SUPPORT_TOL_H = 1.0e-6      # 支持領域の境界 (0.30 H) の float32 の丸めの余裕。境界近傍の φ は (1 − t)³ 程度で丸めより小さい


def _mesher_local_consts():
    """メッシャの局所化の定数 (版・φ = 0 の距離)。取れなければ (None, None)。"""
    try:
        if str(T.DESIGN) not in sys.path:
            sys.path.insert(0, str(T.DESIGN))
        from forge_design.meshing import mesh_sern3d
        return mesh_sern3d.TE_WAKE_MODE_VERSION, float(mesh_sern3d.TE_WAKE_LOCAL_D_ZERO_H)
    except Exception:  # noqa: BLE001
        return None, None


def quality(c, h, sgn=None):
    """各ヘキサの (最小スケール済みヤコビアン [向き sgn]・AR [12 辺の最長/最短]・非正/非有限の頂点ヤコビアンを持つか, sgn)。"""
    J, Js = T.corner_jacobians(c, h)
    if sgn is None:
        sgn = float(np.sign(np.median(J)))
    s = (Js * sgn).min(axis=1)
    bad = (~np.isfinite(J) | ~np.isfinite(Js) | (sgn * J <= 0.0)).any(axis=1)
    del J, Js
    X = c[h]
    L = np.stack([np.linalg.norm(X[:, i] - X[:, j], axis=1) for i, j in EDGES], axis=1)
    del X
    ar = L.max(axis=1) / np.maximum(L.min(axis=1), 1e-300)
    return s, ar, bad, sgn


def support_check(GA, GB, moved):
    """(4 の一部) 局所化の支持領域の外の座標がビット不変か。戻り = dict (status: PASS / FAIL / 判定不能)。"""
    ib = GB["info"]
    r = {"applicable": ib.get("te_wake_mode") == "local"}
    if not r["applicable"]:
        r.update(status=T.UND, note=f"B の te_wake_mode が local でない ({ib.get('te_wake_mode')!r})")
        return r
    ver_now, dzero = _mesher_local_consts()
    r.update(mode_version_B=ib.get("te_wake_mode_version"), mode_version_mesher=ver_now)
    if ver_now is None or ib.get("te_wake_mode_version") != ver_now:
        r.update(status=T.UND, note="B の te_wake_mode_version が現行のメッシャの定数と違う・取れない")
        return r
    try:
        S = T.Struct(ib)
        i0, i1 = int(ib["te_wake_i_first"]), int(ib["te_wake_i_last"])
    except (KeyError, SystemExit) as e:
        r.update(status=T.UND, note=f"B の info から構造・変形区間を作れない ({e})")
        return r
    H = float(GA["H"])
    cA = np.asarray(GA["coords"], float)
    G = S.grid(cA)                                                   # (ni, NJ, nz, 3)
    d = np.abs(G[i0:i1 + 1, :, :, 1] - G[i0:i1 + 1, S.jm:S.jm + 1, :, 1]) / H
    sup = np.zeros(cA.shape[0], bool)
    blk = np.zeros((S.ni, S.NJ, S.nz), bool)
    blk[i0:i1 + 1] = d < dzero + SUPPORT_TOL_H
    sup[:S.N_base] = blk.reshape(-1)
    out = moved & ~sup
    r.update(d_zero_H=dzero, tol_H=SUPPORT_TOL_H, n_support_nodes=int(sup.sum()), n_moved=int(moved.sum()),
             n_moved_outside_support=int(out.sum()), status=T.PASS if not out.any() else T.FAIL)
    if out.any():
        k = int(np.flatnonzero(out)[0])
        r["example_outside"] = {"node": k, "xyz_over_H": [round(float(v) / H, 6) for v in cA[k]]}
    return r


def compare(GA, GB, thr=0.65, ar_max=5000.0):
    """A・B (te_wake_grid_check.load / from_mesher の戻り) の対応するヘキサの比較。戻り = dict (接続・節点数が違えば same_* だけ)。"""
    cA, cB = np.asarray(GA["coords"], float), np.asarray(GB["coords"], float)
    hA, hB = np.asarray(GA["hexes"]), np.asarray(GB["hexes"])
    H = GA["H"]
    out = {"A": GA.get("label"), "B": GB.get("label"), "thr": thr, "ar_max": ar_max}
    out["same_hexes"] = bool(hA.shape == hB.shape and np.array_equal(hA, hB))
    out["same_n_nodes"] = bool(cA.shape == cB.shape)
    if not (out["same_hexes"] and out["same_n_nodes"]):
        return out
    moved = np.any(cA != cB, axis=1)
    reg = moved[hA].any(axis=1)
    sA, arA, badA, sgn = quality(cA, hA)
    sB, arB, badB, _ = quality(cB, hB, sgn)
    new_low = (sA >= thr) & (sB < thr)
    worse_low = (sA < thr) & (sB < sA)
    new_ar = (arA <= ar_max) & (arB > ar_max)
    outside_diff = (~reg) & ((sA != sB) | (arA != arB))
    cen = cB[hB].mean(axis=1) / H

    def loc(mask, key, n=5):
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            return []
        idx = idx[np.argsort(key[idx])][:n]
        return [{"hex": int(k), "A": float(sA[k]), "B": float(sB[k]), "centroid_over_H": [round(float(v), 5) for v in cen[k]]} for k in idx]
    out.update({
        "n_hex": int(hA.shape[0]), "n_moved_nodes": int(moved.sum()), "n_region_hex": int(reg.sum()),
        "min_A": float(sA.min()), "min_B": float(sB.min()),
        "n_lt_thr_A": int((sA < thr).sum()), "n_lt_thr_B": int((sB < thr).sum()),
        "region_n_lt_thr_A": int((reg & (sA < thr)).sum()), "region_n_lt_thr_B": int((reg & (sB < thr)).sum()),
        "region_min_A": float(sA[reg].min()) if reg.any() else None, "region_min_B": float(sB[reg].min()) if reg.any() else None,
        "new_lt_thr": int(new_low.sum()), "worsened_existing_lt_thr": int(worse_low.sum()),
        "outside_region_changed_hex": int(outside_diff.sum()),
        "ar_gt5000_A": int((arA > 5000).sum()), "ar_gt5000_B": int((arB > 5000).sum()), "ar_max_A": float(arA.max()), "ar_max_B": float(arB.max()),
        "new_ar_gt_max": int(new_ar.sum()),
        "region_ar_max_B": float(arB[reg].max()) if reg.any() else None,
        "nonpos_jac_hex_A": int(badA.sum()), "nonpos_jac_hex_B": int(badB.sum()), "orientation_sign_A": sgn,
        "examples_new_lt_thr": loc(new_low, sB), "examples_worsened": loc(worse_low, sB - sA),
    })
    # AR > 5000 の位置 (A)
    idx = np.flatnonzero(arA > 5000)
    if idx.size:
        cx, cy = cen[idx, 0], cen[idx, 1]
        out["ar_gt5000_A_xH_range"] = [float(cx.min()), float(cx.max())]; out["ar_gt5000_A_yH_range"] = [float(cy.min()), float(cy.max())]
        out["ar_gt5000_A_in_region"] = int(reg[idx].sum())
        hist, edges = np.histogram(cx, bins=[-1, 0, 0.5, 1.2, 1.5, 2.2, 5, 12])
        out["ar_gt5000_A_xH_hist"] = list(zip([float(e) for e in edges[:-1]], hist.tolist()))
    out["support"] = support_check(GA, GB, moved)
    return out


def rows_4c(out):
    """plan §5.1 #4c の判定のうちこの道具が出す行。戻り = (総合, [(行, 状態, 値)])。"""
    thr = out["thr"]
    rows = [("(1) 0 格子との対応比較: 新しく 0.65 未満", T.PASS if out["new_lt_thr"] == 0 else T.FAIL, out["new_lt_thr"]),
            ("(1) 既存の 0.65 未満の悪化", T.PASS if out["worsened_existing_lt_thr"] == 0 else T.FAIL, out["worsened_existing_lt_thr"]),
            (f"(1) 新しい AR 超過 (> {out['ar_max']:g})", T.PASS if out["new_ar_gt_max"] == 0 else T.FAIL, out["new_ar_gt_max"])]
    rm = out["region_min_B"]
    rows.append((f"(2) B の変形領域の最小スケール済みヤコビアン ≥ {thr:g}",
                 T.UND if rm is None else (T.PASS if rm >= thr else T.FAIL), rm))
    rows.append(("(2) 全域の非正/非有限の頂点ヤコビアンを持つヘキサ (B)", T.PASS if out["nonpos_jac_hex_B"] == 0 else T.FAIL,
                 out["nonpos_jac_hex_B"]))
    sp = out["support"]
    rows.append(("(4) 支持領域の外の座標のビット不変", sp["status"],
                 sp.get("n_moved_outside_support", sp.get("note"))))
    return T.combine(s for _, s, _ in rows), rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("A"); ap.add_argument("B")
    ap.add_argument("--thr", type=float, default=0.65)
    ap.add_argument("--ar-max", type=float, default=5000.0, help="新しい AR 超過の閾値 (判定 (1))")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    GA, GB = T.load(a.A), T.load(a.B)
    out = compare(GA, GB, a.thr, a.ar_max)
    out["A"], out["B"] = a.A, a.B
    if not (out["same_hexes"] and out["same_n_nodes"]):
        print("UNDECIDABLE: 接続または節点数が違う"); print(json.dumps(out)); return 2
    for k, v in out.items():
        if not isinstance(v, (list, dict)):
            print(f"{k}: {v}")
    print("examples_new_lt_thr:", out["examples_new_lt_thr"][:3]); print("examples_worsened:", out["examples_worsened"][:3])
    print("ar_gt5000_A_xH_hist:", out.get("ar_gt5000_A_xH_hist"))
    print("support:", json.dumps(out["support"], ensure_ascii=False))
    v = "LOCAL_DEGRADATION" if (out["new_lt_thr"] or out["worsened_existing_lt_thr"]) else "NO_LOCAL_DEGRADATION"
    print("HEX_COMPARE VERDICT:", v)
    v4, rows = rows_4c(out)
    out["te_wake_local_rows"] = [{"row": r, "status": s, "value": val} for r, s, val in rows]
    out["te_wake_local_verdict"] = T.VERDICT_WORD[v4]
    print("plan tooling-sern-te-wake-grid §5.1 #4c (この道具の行。(3) と (4) の品質・双対は te_wake_grid_check.py --admission):")
    for r, s, val in rows:
        print(f"  {s:6s} | {r} | {val}")
    print("TE_WAKE_LOCAL VERDICT:", T.VERDICT_WORD[v4])
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
