#!/usr/bin/env python3
r"""温度再現 A/B の判定位置 x=0 での後処理離散化の不確かさ (result レビュー M3、2026-09-29)。

合成場 (放物速度・ρ 一定・源項なし) で、上流の断熱壁と加熱区間の壁 (T_in + step) の角 x=0 を含む問題を
`temp_reproduce.solve_T` で N_r = 16/32/64/128 (軸方向も入れ子で細分、gen_mesh と同じ点列) に解き、
x=0 の内部節点 (判定に使う N_r=32 の全内部節点、壁の隣を含む) で N_r=32 と最細 (128) の差を不確かさとする。
step は実データの角の壁温段差 0.0573 K (`run_0014` の TEMP_REPRODUCE.txt)。
**判定**: max |T_32 − T_128| ≤ 0.005/3 K なら x=0 での判定に使える。細分化で単調に減ることも見る。

    python3 temp_reproduce_x0_uncertainty.py
"""
import numpy as np

import importlib.util
from pathlib import Path

import graetz_common as gc

# graetz_common が case/61 を sys.path の先頭に足すので、同名の case/61/gen_mesh.py を拾わないようにパスで読む
_spec = importlib.util.spec_from_file_location("gen_mesh63", Path(__file__).resolve().parent / "gen_mesh.py")
gen_mesh = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(gen_mesh)
from temp_reproduce import Grid, solve_T, UNC

STEP = 0.0573


def xyz_for(nr):
    up, heat, down = gen_mesh.base_widths()
    f = nr / 32
    def ref(w):
        if f >= 1:
            for _ in range(int(round(np.log2(f)))):
                w = np.repeat(w / 2, 2)
            return w
        return w.reshape(-1, 2).sum(axis=1)
    x = np.concatenate([[-gc.L_UP], -gc.L_UP + np.cumsum(ref(up))])
    x = np.concatenate([x, np.cumsum(ref(heat))[0:]])
    x = np.concatenate([x, gc.L_HEAT + np.cumsum(ref(down))])
    x[np.argmin(np.abs(x))] = 0.0
    r = np.linspace(0, gc.R, nr + 1)
    X, Rr = np.meshgrid(x, r, indexing="ij")
    return np.c_[X.ravel(), Rr.ravel(), np.zeros(X.size)]


res = {}
for nr in (16, 32, 64, 128):
    xyz = xyz_for(nr)
    G = Grid(xyz)
    n = len(xyz); r = xyz[:, 1]
    u = 2 * gc.U_M * np.clip(1 - (r / gc.R) ** 2, 0, None)
    wall = np.where((G.xs >= -1e-12) & (G.xs <= gc.L_HEAT + 1e-9), gc.T_IN + STEP, np.nan)
    T = solve_T(G, np.full(n, gc.RHO), u, np.zeros(n), np.zeros((G.nx, G.ny)), np.full(G.ny, gc.T_IN), wall)
    i0 = int(np.argmin(np.abs(G.xs)))
    res[nr] = (G.r, T[i0])
    print(f"N_r={nr:4d}: {G.nx} 列、x=0 の壁の隣 T−T_in = {T[i0, -2]-gc.T_IN:.5f} K")
# 比較位置は**判定に使う N_r=32 の x=0 の内部節点すべて** (壁の隣 r/R=0.96875 を含む。初版は N_r=16 の半径だけで
# 比べていて、実データの判定に入っている壁の隣の節点を落としていた)
r32 = res[32][0][:-1]
def at(nr, rr_q):
    rr, t = res[nr]
    idx = [int(np.argmin(np.abs(rr - q))) for q in rr_q]
    assert max(abs(rr[i] - q) for i, q in zip(idx, rr_q)) < 1e-12
    return t[idx]
d32 = np.abs(at(32, r32) - at(128, r32))
d64 = np.abs(at(64, r32) - at(128, r32))
j = int(np.argmax(d32))
print(f"  N_r=32 の x=0 内部 {len(r32)} 節点: max |T_32 − T_128| = {d32.max():.2e} K (r/R {r32[j]/gc.R:.5f})、"
      f"同位置 |T_64 − T_128| = {d64.max():.2e} K")
print(f"  壁の隣 r/R=0.96875: T_32 {at(32,[r32[-1]])[0]-gc.T_IN:.5f} / T_64 {at(64,[r32[-1]])[0]-gc.T_IN:.5f} / T_128 {at(128,[r32[-1]])[0]-gc.T_IN:.5f} K (T−T_in)")
ok = d32.max() <= UNC
print(f"VERDICT: {'PASS' if ok else 'FAIL'}  (N_r=32 の不確かさ {d32.max():.2e} K、上限 {UNC:.2e} K; 最細 128 も未収束なら下限の見積もり)")
