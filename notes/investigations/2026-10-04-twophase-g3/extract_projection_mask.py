"""plan condensation-two-phase-default §5.1 #4pjg: G3-b の 1 更新の記録 (FORGE_DIAG_TP_UPDATE の tp_update.h5) から射影集合の節点番号を取り出す。
  python3 extract_projection_mask.py <tp_update.h5> <out.npy> [--buffer update|pre]

射影集合 P = 射影の分岐に入った (RZP_input/kind ≥ 0 ⇔ ρQ0 > 0 かつ ρg > 0) かつ射影 (RZP) が Q1 または Q2 を変えた節点
(pj_ab_judge.py の analyse と同じ定義; 既定は更新のバッファ = 1 更新の中の射影)。
出力: 節点番号 (int64、昇順) の npy。g3a_judge.py --mask に渡す。節点数と座標の要約を表示する (G3-a の h5 と同じメッシュか確かめるため)。
"""
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g3b_judge as g3   # noqa: E402

C_Q2, C_Q1 = 2, 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("out")
    ap.add_argument("--buffer", choices=["update", "pre"], default="update",
                    help="update = 1 更新の中の射影 (既定; #4pj の P_A と同じ)、pre = 組立前処理の射影")
    a = ap.parse_args()
    rec = g3.Rec(a.h5, 1.0)
    ph = "upd" if a.buffer == "update" else "pre"
    kind = rec.s(ph, "RZP_input/kind") if "RZP_input/kind" in rec.idx else None
    if kind is None or not np.isfinite(kind).any():
        raise SystemExit(f"{a.h5}: RZP_input/kind missing or all NaN in the {a.buffer} buffer (h5 written without the #4pj slots?)")
    changed = np.zeros(rec.n, bool)
    for c in (C_Q2, C_Q1):
        r = g3.ba(rec, ph, "RZP", c)
        if r is None:
            raise SystemExit(f"{a.h5}: RZP component {g3.COMPS[c]} not recorded in the {a.buffer} buffer")
        b, af = r
        changed |= (g3.z(af - b) != 0.0)
    P = (kind >= 0) & changed
    idx = np.nonzero(P)[0].astype(np.int64)
    np.save(a.out, idx)
    x = rec.x[idx] if idx.size else np.array([np.nan])
    print(f"{a.h5}: {rec.n} nodes, projection set ({a.buffer} buffer) = {idx.size} nodes "
          f"(changed Q1 {int(np.count_nonzero((kind >= 0) & (g3.z(g3.ba(rec, ph, 'RZP', C_Q1)[1] - g3.ba(rec, ph, 'RZP', C_Q1)[0]) != 0)))}, "
          f"changed Q2 {int(np.count_nonzero((kind >= 0) & (g3.z(g3.ba(rec, ph, 'RZP', C_Q2)[1] - g3.ba(rec, ph, 'RZP', C_Q2)[0]) != 0)))})")
    print(f"  x range of the set [{np.nanmin(x):.6g}, {np.nanmax(x):.6g}] (h5 units); "
          f"node-coordinate checksum sum(ccx) = {float(np.sum(rec.x)):.10e} (compare with g3a_judge's print)")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    sys.exit(main())
