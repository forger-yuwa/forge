#!/usr/bin/env python3
"""G2' : 合併 stencil LSQ の二次場の収束次数 (plan boundary-node-periodic-gradient-fix §6 G2')。

三重周期の立方体 [0,2π]^3 を N = 16/32/64 で作り、節点を ±0.2h の決定論的擬似乱数で動かす (mkmesh.py、周期像は同じ量)。
二次場を継ぎ目中心の局所座標 w で焼き (折返し不連続に触れる節点は除外)、1 step の GPU 勾配 (res_1) を解析勾配と比べる。

  φ (Ux) = w·a + ½κ (w_x² + ½w_y² − 0.7 w_x w_z + 0.3 w_y w_z)        a = (0.3, −0.5, 0.2)
  ψ (Uy) = w·(a_x, a_y, 0) + ½κ (w_x² + ½w_y²)        ∂ψ/∂z ≡ 0 (ゼロ成分)
判定 (§6): e_floor = 10 ε S。継ぎ目・内部の誤差がともに床を超える水準だけで次数を計算し、使える水準が 2 未満なら κ を 10 倍。
継ぎ目誤差の次数 ≥ 0.9、各水準で 継ぎ目誤差/内部誤差 ≤ 2。ゼロ成分は絶対誤差 ≤ e_floor。一様格子 (ジッタ 0) は別行で記録。

使い方: python3 g2p_jitter.py [--kappa 0.2] [--levels 16,32,64]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G   # noqa: E402
import mkmesh          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = "/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad"
A = np.array([0.3, -0.5, 0.2])


def fields(w, kappa):
    x, y, z = w.T
    phi = w @ A + 0.5 * kappa * (x * x + 0.5 * y * y - 0.7 * x * z + 0.3 * y * z)
    gphi = np.stack([A[0] + kappa * (x - 0.35 * z), A[1] + kappa * (0.5 * y + 0.15 * z),
                     A[2] + kappa * (-0.35 * x + 0.15 * y)], axis=1)
    psi = w @ np.array([A[0], A[1], 0.0]) + 0.5 * kappa * (x * x + 0.5 * y * y)
    gpsi = np.stack([A[0] + kappa * x, A[1] + kappa * 0.5 * y, 0 * x], axis=1)
    return phi, gphi, psi, gpsi


def one_level(N, jitter, kappa, tag):
    md = os.path.join(SCR, "g_harness_mesh", f"box{N}_j{jitter:g}")
    h5, bc, q = (os.path.join(md, "mesh.h5"), os.path.join(md, "bcondConfig.yaml"), None)
    if not os.path.exists(h5):
        h5, bc, q = mkmesh.make_box_h5(md, N, jitter)
    else:
        q = [l for l in open(os.path.join(md, "quality.txt")).read().splitlines() if "VERDICT" in l][-1]
    run = os.path.join(SCR, f"g_harness_g2p_{tag}_N{N}_j{jitter:g}")
    cfg = G.base_solver_cfg(sst=False, tracer=False)
    h5r = G.make_run(run, h5, cfg, bc)
    msh = G.Mesh(h5r, os.path.join(run, "bcondConfig.yaml"))
    w, clean, _ = G.wrapped_coords(msh)
    phi, gphi, psi, gpsi = fields(w, kappa)
    G.bake(h5r, G.prim_state_fields(msh, 1.0, phi, psi, 0.0 * phi, 101325.0))
    G.run_forge(run)
    s0 = G.read_res(os.path.join(run, "res_0.h5"), ["Ux", "Uy"])
    g1 = G.read_res(os.path.join(run, "res_1.h5"), [f"dU{v}d{c}" for v in "xy" for c in "xyz"])
    gU = np.stack([g1[f"dUxd{c}"] for c in "xyz"], 1).astype(np.float64)
    gV = np.stack([g1[f"dUyd{c}"] for c in "xyz"], 1).astype(np.float64)
    # 状態が焼いた値のままか (G2' は解析値との比較なので必須)
    same = np.abs(s0["Ux"] - phi.astype(np.float32)) <= 2 * np.spacing(np.abs(phi.astype(np.float32)))
    sel = clean & same
    S = np.linalg.norm(gphi[sel], axis=1).max()
    e = np.abs(gU - gphi).max(axis=1)
    ez = np.abs(gV[:, 2])                                      # ψ のゼロ成分
    seam = msh.nmember > 1
    ref = G.lsq_merged_ref(msh, s0["Ux"])
    g2 = np.abs(gU - ref).max() / np.linalg.norm(ref[sel], axis=1).max()
    return dict(N=N, h=2 * np.pi / N, S=S, e_seam=e[sel & seam].max(), e_int=e[sel & ~seam].max(),
                ez_seam=ez[sel & seam].max(), ez_int=ez[sel & ~seam].max(), n_seam=int((sel & seam).sum()),
                n_int=int((sel & ~seam).sum()), g2=g2, run=run, q=q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=0.2)
    ap.add_argument("--levels", default="16,32,64")
    a = ap.parse_args()
    levels = [int(x) for x in a.levels.split(",")]
    out = ["# G2': 合併 stencil LSQ の二次場の収束 (plan boundary-node-periodic-gradient-fix §6 G2')",
           f"harness revision: {G.git_rev()}", ""]
    kappa = a.kappa
    for attempt in range(2):
        rows = [one_level(N, 0.2, kappa, f"k{kappa:g}") for N in levels]
        ef = [10 * G.EPS32 * r["S"] for r in rows]
        usable = [i for i, r in enumerate(rows) if r["e_seam"] > ef[i] and r["e_int"] > ef[i]]
        if len(usable) >= 2:
            break
        out.append(f"(κ={kappa:g}: 床を超える水準が {len(usable)} → κ を 10 倍して再試験)")
        kappa *= 10
    out.append(f"ジッタ ±0.2h、κ = {kappa:g}、a = {tuple(A)}。e = max 成分誤差 (対 解析勾配)、e_floor = 10 ε S")
    out.append("| N | h | 継ぎ目節点 | 内部節点 | S | e_floor | 継ぎ目 e | 内部 e | 継ぎ目/内部 | 次数 (継ぎ目) | 次数 (内部) | ゼロ成分 継ぎ目 / 内部 | G2 (GPU−CPU double)/S | run |")
    out.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    ok = True
    prev = None
    orders = []
    for i, r in enumerate(rows):
        ratio = r["e_seam"] / r["e_int"]
        ok &= ratio <= 2.0
        o_s = o_i = "-"
        if prev is not None and i in usable and prev in usable:
            p_s = np.log(rows[prev]["e_seam"] / r["e_seam"]) / np.log(rows[prev]["h"] / r["h"])
            p_i = np.log(rows[prev]["e_int"] / r["e_int"]) / np.log(rows[prev]["h"] / r["h"])
            orders.append(p_s); o_s, o_i = f"{p_s:.2f}", f"{p_i:.2f}"
        prev = i
        zok = r["ez_seam"] <= ef[i] and r["ez_int"] <= ef[i]
        out.append(f"| {r['N']} | {r['h']:.4f} | {r['n_seam']} | {r['n_int']} | {r['S']:.3f} | {ef[i]:.1e} | {r['e_seam']:.3e} | "
                   f"{r['e_int']:.3e} | {ratio:.2f} | {o_s} | {o_i} | {r['ez_seam']:.1e} / {r['ez_int']:.1e} {'ok' if zok else 'NG'} | "
                   f"{r['g2']:.1e} | {os.path.basename(r['run'])} |")
    ok_order = len(orders) >= 1 and min(orders) >= 0.9
    zero_ok = all(r["ez_seam"] <= ef[i] and r["ez_int"] <= ef[i] for i, r in enumerate(rows))
    out.append(f"メッシュ品質: " + "; ".join(f"N={r['N']}: {r['q']}" for r in rows))
    # 一様格子 (ジッタ 0) の別行
    u = one_level(32, 0.0, kappa, f"k{kappa:g}")
    efu = 10 * G.EPS32 * u["S"]
    out.append(f"\n一様格子 (ジッタ 0、N=32): 継ぎ目 e = {u['e_seam']:.2e}、内部 e = {u['e_int']:.2e}、e_floor = {efu:.1e} → "
               f"{'床以下 (対称 stencil で二次場の勾配は厳密)' if max(u['e_seam'], u['e_int']) <= efu else '床を超える'}。"
               f"ゼロ成分 {max(u['ez_seam'], u['ez_int']):.1e}。run {os.path.basename(u['run'])}")
    out.append(f"\n次数 (継ぎ目) ≥ 0.9: {'ok' if ok_order else 'NG'} ({', '.join(f'{p:.2f}' for p in orders)})。"
               f"継ぎ目/内部 ≤ 2: {'ok' if ok else 'NG'}。ゼロ成分 ≤ e_floor: {'ok' if zero_ok else 'NG'}")
    out.append(f"VERDICT G2' (次数・比): {'PASS' if (ok_order and ok) else 'FAIL'}")
    out.append(f"VERDICT G2' (ゼロ成分 ≤ e_floor): {'PASS' if zero_ok else 'FAIL'}")
    txt = "\n".join(out) + "\n"
    p = os.path.join(HERE, "G2p_jitter.txt")
    open(p, "w").write(txt)
    print(txt)
    print("->", p)


if __name__ == "__main__":
    main()
