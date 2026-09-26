#!/usr/bin/env python3
"""G2' : 合併 stencil LSQ の二次場の収束次数 (plan boundary-node-periodic-gradient-fix §6 G2'、§5.1 #6b)。

三重周期の立方体 [0,2π]^3 を N = 16/32/64/128 で作り、節点を ±0.2h の決定論的擬似乱数で動かす (mkmesh.py、周期像は同じ量)。
二次場を継ぎ目中心の局所座標 w で焼き (折返し不連続に触れる節点は除外)、1 step の GPU 勾配 (res_1) を解析勾配と比べる。

  φ (Ux) = w·a + ½κ (w_x² + ½w_y² − 0.7 w_x w_z + 0.3 w_y w_z)        a = (0.3, −0.5, 0.2)
  ψ (Uy) = w·(a_x, a_y, 0) + ½κ (w_x² + ½w_y²)        ∂ψ/∂z ≡ 0 (ゼロ成分)

判定 (§6 G2' の 2026-09-26 訂正と §5.1 #6b の「測る前に固定」):
  - 誤差は成分ごと (φ_x, φ_y, φ_z, ψ_x, ψ_y, ψ_z = ゼロ成分) の max |∇_gpu − ∇_解析| を継ぎ目 (member≥2)・内部で取る。
  - e_floor = 10 ε S。継ぎ目・内部の誤差がともに床を超える水準だけで次数を計算し、使える水準が 2 未満なら κ を 10 倍。
  - **ジッタ格子ではゼロ成分も他成分と同じ判定** (次数 ≥ 0.9・継ぎ目/内部 ≤ 2)。「ゼロ成分 ≤ e_floor」は一様格子の行のみ。
  - 本判定: 最細対 64→128 の継ぎ目の次数が**全成分**で ≥ 0.9、かつ 128³ で継ぎ目/内部 ≤ 2 (全成分)。
    32→64 は補助 (表示)、16→32 は漸近域外として記録のみ。64→128 が 0.9 未満なら「判定不能」とし、内部の次数を並べる。
  - 一様格子 (ジッタ 0、N=32) は別行: 床以下 = 厳密、ゼロ成分 ≤ e_floor (記録)。

使い方: python3 g2p_jitter.py [--kappa 0.2] [--levels 16,32,64,128] [--scratch DIR]  (scratch は環境変数 G_HARNESS_SCRATCH でも可)
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G   # noqa: E402
import mkmesh          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.environ.get("G_HARNESS_SCRATCH", "/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad")
COMPS = ["φ_x", "φ_y", "φ_z", "ψ_x", "ψ_y", "ψ_z (ゼロ)"]
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
    E = np.concatenate([np.abs(gU - gphi), np.abs(gV - gpsi)], axis=1)       # (n, 6) 成分ごとの誤差
    seam = msh.nmember > 1
    ref = G.lsq_merged_ref(msh, s0["Ux"])
    g2 = np.abs(gU - ref).max() / np.linalg.norm(ref[sel], axis=1).max()
    return dict(N=N, h=2 * np.pi / N, S=S, es=E[sel & seam].max(axis=0), ei=E[sel & ~seam].max(axis=0),
                n_seam=int((sel & seam).sum()), n_int=int((sel & ~seam).sum()), g2=g2, run=run, q=q)


def order(r0, r1, key, c):
    return np.log(r0[key][c] / r1[key][c]) / np.log(r0["h"] / r1["h"])


def main():
    global SCR
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=0.2)
    ap.add_argument("--levels", default="16,32,64,128")
    ap.add_argument("--scratch", default=SCR)
    a = ap.parse_args()
    SCR = a.scratch
    levels = [int(x) for x in a.levels.split(",")]
    out = ["# G2': 合併 stencil LSQ の二次場の収束 (plan boundary-node-periodic-gradient-fix §6 G2'、§5.1 #6b)",
           f"harness revision: {G.git_rev()}", ""]
    kappa = a.kappa
    for attempt in range(2):
        rows = [one_level(N, 0.2, kappa, f"k{kappa:g}") for N in levels]
        ef = [10 * G.EPS32 * r["S"] for r in rows]
        # 使える水準 = 全成分で継ぎ目・内部の誤差がともに床を超える水準
        usable = [i for i, r in enumerate(rows) if (r["es"] > ef[i]).all() and (r["ei"] > ef[i]).all()]
        if len(usable) >= 2:
            break
        out.append(f"(κ={kappa:g}: 床を超える水準が {len(usable)} → κ を 10 倍して再試験)")
        kappa *= 10
    out.append(f"ジッタ ±0.2h、κ = {kappa:g}、a = {tuple(float(x) for x in A)}。誤差 = 成分ごとの max |∇_gpu − ∇_解析|、e_floor = 10 ε S")
    out.append("")
    out.append("## 水準ごと (継ぎ目 / 内部 / 比)")
    out.append("| N | h | 継ぎ目節点 | 内部節点 | S | e_floor | " + " | ".join(f"{c} 継ぎ目 / 内部 (比)" for c in COMPS)
               + " | G2 (GPU−CPU double)/S | run |")
    out.append("| --- " * (8 + len(COMPS)) + "|")
    for i, r in enumerate(rows):
        cells = [f"{r['es'][c]:.3e} / {r['ei'][c]:.3e} ({r['es'][c] / r['ei'][c]:.2f})" for c in range(len(COMPS))]
        out.append(f"| {r['N']} | {r['h']:.4f} | {r['n_seam']} | {r['n_int']} | {r['S']:.3f} | {ef[i]:.1e} | "
                   + " | ".join(cells) + f" | {r['g2']:.1e} | {os.path.basename(r['run'])} |")
    out.append(f"メッシュ品質: " + "; ".join(f"N={r['N']}: {r['q']}" for r in rows))
    out.append("")
    out.append("## 次数 (継ぎ目 / 内部)。16→32 は記録のみ、32→64 は補助、64→128 が本判定")
    out.append("| 対 | 扱い | " + " | ".join(COMPS) + " | 継ぎ目 最小 |")
    out.append("| --- " * (3 + len(COMPS)) + "|")
    pair_orders = {}
    for i in range(1, len(rows)):
        r0, r1 = rows[i - 1], rows[i]
        tag = f"{r0['N']}→{r1['N']}"
        role = {16: "記録のみ", 32: "補助", 64: "本判定"}.get(r0["N"], "-")
        if (i - 1) in usable and i in usable:
            ps = [order(r0, r1, "es", c) for c in range(len(COMPS))]
            pi = [order(r0, r1, "ei", c) for c in range(len(COMPS))]
            pair_orders[(r0["N"], r1["N"])] = (ps, pi)
            out.append(f"| {tag} | {role} | " + " | ".join(f"{ps[c]:.3f} / {pi[c]:.3f}" for c in range(len(COMPS)))
                       + f" | {min(ps):.3f} |")
        else:
            out.append(f"| {tag} | {role} | " + " | ".join("床以下で計算しない" for _ in COMPS) + " | - |")
    # 本判定
    last = rows[-1]
    fine = pair_orders.get((64, 128))
    ratio128 = last["es"] / last["ei"] if last["N"] == 128 else None
    out.append("")
    if fine is None or ratio128 is None:
        verdict = "判定不能 (64→128 の対が無い、または床以下)"
    else:
        ps, pi = fine
        rok = bool((ratio128 <= 2.0).all())
        if min(ps) < 0.9:
            verdict = ("判定不能 (64→128 の継ぎ目の次数 最小 %.3f < 0.9)。内部の次数: %s。128³ の継ぎ目/内部: %s"
                       % (min(ps), ", ".join(f"{COMPS[c]} {pi[c]:.3f}" for c in range(len(COMPS))),
                          ", ".join(f"{COMPS[c]} {ratio128[c]:.2f}" for c in range(len(COMPS)))))
        else:
            verdict = ("PASS" if rok else "FAIL") + (" (64→128 の継ぎ目の次数 最小 %.3f ≥ 0.9、128³ の継ぎ目/内部 最大 %.2f %s 2)"
                                                     % (min(ps), ratio128.max(), "≤" if rok else ">"))
    aux = pair_orders.get((32, 64))
    if aux is not None:
        out.append(f"補助 32→64: 継ぎ目の次数 最小 {min(aux[0]):.3f} ({'≥' if min(aux[0]) >= 0.9 else '<'} 0.9)。"
                   f"各水準の継ぎ目/内部 最大: " + ", ".join(f"N={r['N']} {(r['es'] / r['ei']).max():.2f}" for r in rows))
    # 一様格子 (ジッタ 0) の別行
    u = one_level(32, 0.0, kappa, f"k{kappa:g}")
    efu = 10 * G.EPS32 * u["S"]
    ez = max(u["es"][5], u["ei"][5])
    out.append(f"\n一様格子 (ジッタ 0、N=32、記録): 継ぎ目 e = {u['es'].max():.2e}、内部 e = {u['ei'].max():.2e}、e_floor = {efu:.1e} → "
               f"{'床以下 (対称 stencil で二次場の勾配は厳密)' if max(u['es'].max(), u['ei'].max()) <= efu else '床を超える'}。"
               f"ゼロ成分 {ez:.1e} ({'≤' if ez <= efu else '>'} e_floor)。run {os.path.basename(u['run'])}")
    out.append(f"\nVERDICT G2' (ジッタ、64→128 の継ぎ目の次数 全成分 ≥ 0.9 かつ 128³ の継ぎ目/内部 ≤ 2): {verdict}")
    out.append(f"VERDICT G2' (一様格子 N=32、ゼロ成分 ≤ e_floor、記録): {'PASS' if ez <= efu else 'FAIL'}")
    txt = "\n".join(out) + "\n"
    p = os.path.join(HERE, "G2p_jitter.txt")
    open(p, "w").write(txt)
    print(txt)
    print("->", p)


if __name__ == "__main__":
    main()
