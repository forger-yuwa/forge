"""plan architecture-float-state-double-geometry §6.57 (事前登録、codex diagnose `2026-10-12-speedup-tests-diagnose.md`) の判定:
前処理の経路への切り替えが node の境界の更新を壊さないかの 1 step の確認。G0・Q32 の起点・point・FP64 (`~/forge-fgeom11-fp64`)、両腕とも
`precondEps 1.0` (β = 1 で前処理は数式上は無効)、変えるのは `lowMachPrecond` 0 (P0) / 3 (P3) だけ、各 2 本・1 step、`FORGE_DIAG_STAGE_DUMP`。
量: stage.h5 の d (commit が足す要求の更新、境界条件をかけ直す前)、R (solve の直前の右辺)、Q0。成分 ro・roUx・roUy・roe (roUz は 2D で 0 なら除く)。
領域: 軸 (axis_flag)、壁 (wall_flag。等温壁は iso_wall_flag で別記)、その他の境界 (bnode_flag で軸・壁でない = 入口・出口)、内部 (どれでもない)。
指標: 相対差 e = ‖d_P3 − d_P0‖ / ‖d_P0‖ (各腕の 1 本目どうし)、再実行の幅 s = max(‖d_P0a − d_P0b‖, ‖d_P3a − d_P3b‖) / ‖d_P0‖、大きさの比 g = ‖d_P3‖ / ‖d_P0‖。
前提 (欠ければ判別不能): 4 本の Q0 がビット一致、R が腕の間で再実行の幅の中 (RHS は変えない設定)、stage.h5 がそろい path の属性が 4 本とも
  blockDPLUR (ライン無し)、commit の自己照合 (q == fl(b + d)) の不一致 0、各領域で ‖d_P0‖ > 0 (0 の成分は拘束として別に見る)。
結果 A (切り替えは純粋な前処理の操作ではない): ある領域・成分で e > max(10 s, 0.1)、または g ≤ 0.01、または拘束の成分 (P0 で d が全節点 0) が P3 で 0 でない。
結果 B: すべての領域・成分で e ≤ max(s, 0.01) かつ拘束の成分の違反なし。それ以外は判別不能。→ _band_ab/cold_pair/rlg_judge.json"""
import json
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
P0 = ("run_0756_lg_p0a", "run_0757_lg_p0b"); P3 = ("run_0758_lg_p3a", "run_0759_lg_p3b")
Q = ("ro", "roUx", "roUy", "roUz", "roe")
OUT = {"undecidable": [], "pre": {}, "table": {}, "constraints": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def load(r):
    p = HERE / r / "stage.h5"
    if not p.exists(): return None
    with h5py.File(p, "r") as h:
        d = {g: {q: np.asarray(h[f"/{g}/{q}"][:], np.float64) for q in Q} for g in ("Q0", "R", "b", "d", "q")}
        fl = {k: np.asarray(h[f"/cells/{k}_flag"][:]) for k in ("bnode", "axis", "wall", "iso_wall") if f"/cells/{k}_flag" in h}
        at = {k: (h.attrs[k].decode() if isinstance(h.attrs[k], bytes) else (np.asarray(h.attrs[k]).item() if np.asarray(h.attrs[k]).size == 1 else str(h.attrs[k]))) for k in h.attrs}
    return d, fl, at
D = {r: load(r) for r in P0 + P3}
for r, v in D.items():
    if v is None: und(f"{r}: stage.h5 が無い (経路が診断の対象外の可能性)")
if not OUT["undecidable"]:
    for r, (d, fl, at) in D.items():
        OUT["pre"][f"{r} 属性"] = {k: str(at[k])[:80] for k in ("path", "commit_kernel", "d_source", "step") if k in at}
        if "blockDPLUR" not in str(at.get("path", "")) or "line" in str(at.get("path", "")).lower(): und(f"{r}: path が blockDPLUR (ライン無し) でない ({at.get('path')})")
        mism = [q for q in Q if not np.array_equal(d["q"][q], (d["b"][q] + d["d"][q]))]
        if mism: OUT["pre"][f"{r} q != b + d の量 (double で照合、記録)"] = mism
    d0 = D[P0[0]][0]
    for r in P0[1:] + P3:
        if any(not np.array_equal(D[r][0]["Q0"][q], d0["Q0"][q]) for q in Q): und(f"{r}: Q0 が {P0[0]} とビット一致しない")
if not OUT["undecidable"]:
    fl = D[P0[0]][1]; n = len(D[P0[0]][0]["d"]["ro"])
    ax = fl.get("axis", np.zeros(n)).astype(bool); wa = fl.get("wall", np.zeros(n)).astype(bool); bn = fl.get("bnode", np.zeros(n)).astype(bool)
    REG = {"軸": ax, "壁": wa & ~ax, "その他の境界": bn & ~ax & ~wa, "内部": ~bn & ~ax & ~wa}
    OUT["pre"]["領域の節点数"] = {k: int(v.sum()) for k, v in REG.items()}
    comps = [q for q in Q if not (q == "roUz" and all(np.all(D[r][0]["d"][q] == 0) for r in D))]
    nrm = lambda a: float(np.sqrt(np.sum(a * a)))  # noqa: E731
    # 前提: R が腕の間で再実行の幅の中
    for nm, m in REG.items():
        for q in comps:
            r0a, r0b, r3a = (D[x][0]["R"][q][m] for x in (P0[0], P0[1], P3[0]))
            base = nrm(r0a) or 1e-300
            if nrm(r3a - r0a) / base > max(10 * nrm(r0b - r0a) / base, 1e-6): und(f"R が腕で違う ({nm}・{q}: {nrm(r3a - r0a) / base:.2e})")
    ra, rb = [], []
    for nm, m in REG.items():
        for q in comps:
            a0, b0, a3, b3 = (D[x][0]["d"][q][m] for x in (P0[0], P0[1], P3[0], P3[1]))
            n0 = nrm(a0)
            if n0 == 0:
                viol = float(max(np.max(np.abs(a3)) if a3.size else 0.0, np.max(np.abs(b3)) if b3.size else 0.0))
                OUT["constraints"][f"{nm}|{q}"] = {"P0 は全節点 0": True, "P3 の max|d|": viol}
                ra.append(viol > 0); rb.append(viol == 0); continue
            e = nrm(a3 - a0) / n0; s = max(nrm(a0 - b0), nrm(a3 - b3)) / n0; g = nrm(a3) / n0
            OUT["table"][f"{nm}|{q}"] = dict(e=e, s=s, g=g, n0=n0)
            print(f"{nm:6s} {q:5s} e {e:.3e}  s {s:.3e}  g {g:.4f}  ‖d_P0‖ {n0:.3e}")
            ra.append(e > max(10 * s, 0.1) or g <= 0.01); rb.append(e <= max(s, 0.01))
    for k, v in OUT["constraints"].items(): print("拘束", k, v)
    verdict = ("結果 A: 前処理の経路への切り替えだけで境界・軸の更新が変わる (純粋な前処理の操作ではない。境界の経路を先に直す)" if any(ra) else
               "結果 B: 切り替えだけでは更新はほぼ変わらない (今回の入力で重大な境界の更新の障害は出ない)" if all(rb) else "判別不能")
else:
    verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.57: {verdict}")
(HERE / "_band_ab/cold_pair/rlg_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rlg_judge.json")
