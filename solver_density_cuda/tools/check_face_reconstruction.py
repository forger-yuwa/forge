#!/usr/bin/env python3
"""再構成した面状態を後処理で検査する (node)。**有界性 G1** と **物理性 G3** を分けて出す。

SU2 は MUSCL 再構成の後に `P<0 || rho<0 || 音速^2<0` を検査し、該当エッジを 20 反復 1 次に落とす
(`CEulerSolver.cpp` の `bad_recon`)。**forge にはこの検査が無い**ので、壊れた面状態がそのまま流束へ入る。

## 使い方が難しい点が 3 つある (どれも実際に踏んだ)

**(1) 座標**: `CELLS/centCoords` は**実行時にノード座標へ置換される**
(`main.cpp` の `nodeValueAtNode` → `mesh.cpp`。`forge_run.log` に `centCoords <- node coords ... max centroid shift`)。
本ツールは `res_*.h5` の `MESH/COORD` を使う。メッシュ h5 の `centCoords` を使うと別の点を評価したことになる。

**(2) 時相**: 定常陰解法は「勾配・リミッタを q^n から作る → 流束 → q^{n+1} に更新 → 出力」の順なので、
`res_n.h5` の **保存量は更新後 (q^{n+1})、勾配・リミッタは更新前 (q^n 由来)** である
(`updateVariablesOuter` は勾配を作り直さない)。**同じ snapshot の中で組み合わせてはいけない。**
既定 `--state-from prev` は状態を `res_{n-1}.h5` から、勾配・リミッタを `res_n.h5` から取る
(= 流束評価時の組み合わせ)。`--state-from same` は旧挙動で、**診断目的以外では使わない**。

**(3) 増分の形**: `convMethod: 2` は勾配射影に加えて隣接値差の項を持つ (`interp_MUSCL_3rd`)。
`--conv-method` で指定する (既定 1)。

## 出す判定

- **G1 有界性**: 再構成した面値が**そのノードの近傍全体の min/max** を外れた面側の数 (丸め許容 `--tol`)。
  リミッタが保証するのはこの集合であって、面の 2 セルの範囲ではない (面の 2 セルで測ると Barth でも大量に逸脱に見える)。
  Barth は厳密有界を要求できる。Venkatakrishnan は平滑化項があるので逸脱しうる (逸脱量の比較に使う)。
- **G3 物理性**: 再構成した `P`・`rho` が**有限かつ正**か。SU2 の `bad_recon` と同じ検査。

**入力に非有限値がある / 指定 snapshot が無い / 検査件数が 0 のときは不合格 (終了コード 1)**。
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_lsq_gradient import parse_plane_cells  # noqa: E402

VARS = ("ro", "Ux", "Uy", "Uz", "P")
GRAD = {"ro": ("drodx", "drody", "drodz"), "Ux": ("dUxdx", "dUxdy", "dUxdz"),
        "Uy": ("dUydx", "dUydy", "dUydz"), "Uz": ("dUzdx", "dUzdy", "dUzdz"),
        "P": ("dPdx", "dPdy", "dPdz")}
LIM = {v: f"limiter_{v}" for v in VARS}


def _read(path, keys):
    with h5py.File(path, "r") as f:
        V = f["VALUE"]
        missing = [k for k in keys if k not in V]
        if missing:
            raise SystemExit(f"{path}: 変数不足 {missing}  → `output: {{level: 2}}` で出し直す")
        out = {k: np.asarray(V[k][()]).ravel().astype(np.float64) for k in keys}
        out["_xyz"] = np.asarray(f["MESH/COORD"][()]).reshape(-1, 3).astype(np.float64)
    return out


def scan(state_path, grad_path, i0, i1, conv_method, tol):
    """戻り値: dict(n_faces, g1[var]=逸脱数, g3=非物理数, pmin, nonfinite, worst)"""
    st = _read(state_path, list(VARS))
    gr = _read(grad_path, [g for v in VARS for g in GRAD[v]] + [LIM[v] for v in VARS])
    xyz = st["_xyz"]
    nonfinite = {}
    for k in VARS:
        nb = int((~np.isfinite(st[k])).sum())
        if nb:
            nonfinite[k] = nb
    for v in VARS:
        for g in GRAD[v] + (LIM[v],):
            nb = int((~np.isfinite(gr[g])).sum())
            if nb:
                nonfinite[g] = nb
    d0 = 0.5 * (xyz[i1] - xyz[i0])
    res = {"g1": {}, "g3": 0, "pmin": np.inf, "nonfinite": nonfinite, "worst": None,
           "n_faces": len(i0)}
    if nonfinite:
        return res
    kk = 1.0 / 3.0
    for v in VARS:
        q = st[v]
        gx, gy, gz = (gr[g] for g in GRAD[v])
        psi = gr[LIM[v]]

        def rec(i, dv, other):
            d = psi[i] * (gx[i] * dv[:, 0] + gy[i] * dv[:, 1] + gz[i] * dv[:, 2])
            if conv_method == 2:
                d = psi[i] * (0.5 * kk * (q[other] - q[i])) + (1.0 - kk) * d
            return q[i] + d

        L, R = rec(i0, d0, i1), rec(i1, -d0, i0)
        # リミッタが保証するのは「その**ノードの近傍全体**の min/max」であって面の 2 セルの範囲ではない
        # (limiter_d.cu の pass1 が近傍 min/max を作る)。同じ集合で判定しないと Barth でも大量に逸脱に見える。
        nmax = q.copy(); nmin = q.copy()
        np.maximum.at(nmax, i0, q[i1]); np.maximum.at(nmax, i1, q[i0])
        np.minimum.at(nmin, i0, q[i1]); np.minimum.at(nmin, i1, q[i0])
        # 丸め許容の尺度は「近傍レンジ」と「値の大きさ」の**大きい方**。
        # レンジだけだと自由流 (レンジがちょうど 0) で発散し、値だけだと 0 を跨ぐ速度で許容が潰れる。
        magL = np.maximum(np.abs(nmax[i0]), np.abs(nmin[i0]))
        magR = np.maximum(np.abs(nmax[i1]), np.abs(nmin[i1]))
        scL = np.maximum(nmax[i0] - nmin[i0], magL) * tol + 1e-300
        scR = np.maximum(nmax[i1] - nmin[i1], magR) * tol + 1e-300
        res["g1"][v] = int(((L > nmax[i0] + scL) | (L < nmin[i0] - scL)
                            | (R > nmax[i1] + scR) | (R < nmin[i1] - scR)).sum())
        if v in ("P", "ro"):
            bad = ~np.isfinite(L) | ~np.isfinite(R) | (L <= 0) | (R <= 0)
            res["g3"] += int(bad.sum())
            m = np.nanmin(np.minimum(L, R))
            if v == "P":
                res["pmin"] = m
                k = int(np.argmin(np.where(np.isfinite(np.minimum(L, R)), np.minimum(L, R), np.inf)))
                res["worst"] = (int(i0[k]), int(i1[k]), xyz[i0[k]])
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="再構成した面状態の有界性 (G1) と物理性 (G3) を検査する")
    ap.add_argument("run_dir")
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--steps", default=None, help="カンマ区切り (既定: 全 res_*.h5)")
    ap.add_argument("--conv-method", type=int, default=1, choices=[0, 1, 2])
    ap.add_argument("--state-from", default="prev", choices=["prev", "same"],
                    help="保存量をどの snapshot から取るか。既定 prev = 流束評価時の組み合わせ (docstring 参照)")
    ap.add_argument("--tol", type=float, default=1e-5,
                    help="有界性の丸め許容 (値の大きさに対する相対)。float32 の丸めを逸脱に数えないため既定 1e-5")
    a = ap.parse_args(argv)
    rd = a.run_dir
    mesh = a.mesh or next(p for p in sorted(glob.glob(f"{rd}/*.h5")) if "/res_" not in p)
    with h5py.File(mesh, "r") as f:
        st = np.asarray(f["PLANES/STRUCT"][()])
        nfaces = np.asarray(f["PLANES/surfArea"][()]).shape[0]
    ic0, ic1 = parse_plane_cells(st, nfaces)
    sel = ic1 >= 0
    i0, i1 = ic0[sel], ic1[sel]
    have = sorted(int(os.path.basename(p)[4:-3]) for p in glob.glob(f"{rd}/res_[0-9]*.h5"))
    steps = [int(x) for x in a.steps.split(",")] if a.steps else have
    print(f"{rd}: 内部面 {int(sel.sum())} / 全面 {nfaces}  conv_method={a.conv_method}  state_from={a.state_from}")
    print(f"  {'step':>6s} {'G1 逸脱':>8s} {'G3 非物理':>9s} {'再構成 P の最小':>16s}   最悪面")
    fail, checked = [], 0
    for s in steps:
        cur = f"{rd}/res_{s}.h5"
        prev = f"{rd}/res_{max([h for h in have if h < s], default=-1)}.h5" if a.state_from == "prev" else cur
        if not os.path.exists(cur):
            fail.append(f"step {s}: res が無い"); continue
        if a.state_from == "prev" and not os.path.exists(prev):
            fail.append(f"step {s}: 直前の res が無い (state-from prev)"); continue
        r = scan(prev, cur, i0, i1, a.conv_method, a.tol)
        checked += 1
        if r["nonfinite"]:
            print(f"  {s:6d} {'-':>8s} {'-':>9s} {'-':>16s}   入力に非有限値 {r['nonfinite']}")
            fail.append(f"step {s}: 入力に非有限値 {r['nonfinite']}")
            continue
        g1 = sum(r["g1"].values())
        w = r["worst"]
        loc = f"{w[0]}->{w[1]} at ({w[2][0]:.5f}, {w[2][1]:.5f}, {w[2][2]:.5f})" if w else ""
        print(f"  {s:6d} {g1:8d} {r['g3']:9d} {r['pmin']:16.2f}   {loc}")
        if r["g3"]:
            fail.append(f"step {s}: 非物理な面状態 {r['g3']}")
    print()
    if checked == 0:
        print("VERDICT: FAIL (検査した snapshot が 0)")
        return 1
    if fail:
        print(f"VERDICT: FAIL ({len(fail)} 件)")
        for m in fail[:8]:
            print(f"  - {m}")
        return 1
    print(f"VERDICT: OK (G3 非物理ゼロ。G1 の逸脱数は上表。検査 {checked} snapshot)")
    print("注意: G1 の逸脱は Barth では 0 を要求できるが Venkatakrishnan では平滑化項の分だけ残る (比較に使う)。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
