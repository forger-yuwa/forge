#!/usr/bin/env python3
"""G2′ 固定形状縮小試験 (plan boundary-node-periodic-gradient-fix §5.1 #8b、codex result-2 M2)。

32³ ジッタ格子 (±0.2h、`g_harness_mesh/box32_j0.2`) の**実在の継ぎ目 group** (面 2・辺 4・角 8 member) について、
group の全 member の incidence (変位 d = 格納 float32 座標の差を double で) を取り出し、無次元形状を固定したまま
d → s·d (s = 1, 1/2, 1/4、2 の冪なので縮小は厳密) と縮小して、二次場 φ = a·x + ½ xᵀA x (原点 = group の root) の
勾配を G2 の参照 `gharness.lsq_merged_ref` (合併 stencil LSQ の double 参照、GPU と ≤ 1e-5·S で一致済み) を
**そのまま**呼んで求める。解析勾配は原点で a なので、誤差 = |∇_ref − a| (成分ごと)。

- 部分メッシュ: `gharness.Mesh` のインスタンスを __new__ で作り、group の member とその隣接だけを持たせる
  (inc_m / inc_j / inc_d / root / nmember)。α (同値類の重複数) は `Mesh.alpha` がそのまま決める (形状不変)。
- 自己検査: s = 1 で実在の格納場 (二次場を実座標で焼いたもの) を与えた部分メッシュの結果が、全メッシュでの
  `lsq_merged_ref` の同じ root の値と一致すること (抽出が stencil を落としていないこと)。
- 内部: member 1 の節点を同じ手順で縮小し、誤差/h の定数 C を並べる。

判定 (§5.1 #8b、測る前に固定): 面 2・辺 4・角 8 member の実在 group 各 ≥ 2 個・計 ≥ 6 個で、h→h/2 と h/2→h/4 の
継ぎ目の次数が**全成分** ≥ 0.9。
**注記**: 形状を固定した縮小では、LSQ の二次場の誤差は s に厳密に比例する (M は s に不変、右辺は s に比例) ので、
次数はほぼ 1 になる。この試験が示すのは合併作用素の 1 次整合 (打ち切り誤差の主要項が h に比例すること) であり、
ジッタ系列 (評価点集合と stencil 形状が格子ごとに変わる、G2p_jitter.txt の 0.855) の判定不能を覆すものではない。

使い方: python3 g2p_fixedshape.py [--scratch DIR]   (メッシュは DIR/g_harness_mesh/box32_j0.2)
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.environ.get("G_HARNESS_SCRATCH", "/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad")
A_LIN = np.array([0.3, -0.5, 0.2])
KAPPA = 0.2
# ½ xᵀ H x = ½κ (x² + ½y² − 0.7xz + 0.3yz) (g2p_jitter.py の φ と同じ二次形)
HESS = KAPPA * np.array([[1.0, 0.0, -0.35], [0.0, 0.5, 0.15], [-0.35, 0.15, 0.0]])
SCALES = (1.0, 0.5, 0.25)
COMPS = ("x", "y", "z")
N_PER_KIND = 3        # 各種別から選ぶ group 数 (等間隔、決定論的)
N_INTERIOR = 6
DFLOOR = 1e-14        # 倍精度丸めの目安 (|a| ~ 0.5 に対し)。全縮小率でこれ以下の成分は次数を定義できないので別記する


def sub_mesh(msh, members):
    """members (global index、先頭が root) の incidence だけを持つ部分メッシュ。戻り値 (sub, P_local, root_local)。
    P_local: 局所座標 (root = 原点)。隣接の位置は member からの変位 inc_d で決める (同じ隣接が複数 member から
    見えるときは 1e-9·|d| で一致を検査)。"""
    mem = list(members)
    es = np.nonzero(np.isin(msh.inc_m, mem))[0]
    nodes = list(mem)
    for j in msh.inc_j[es]:
        if j not in nodes:
            nodes.append(int(j))
    loc = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)
    P = np.full((n, 3), np.nan)
    P[:len(mem)] = 0.0                                   # member は全員 root と同じ点 (周期像)
    for e in es:
        m, j = loc[int(msh.inc_m[e])], loc[int(msh.inc_j[e])]
        pj = P[m] + msh.inc_d[e]
        if np.isnan(P[j, 0]):
            P[j] = pj
        elif np.abs(P[j] - pj).max() > 1e-9 * np.linalg.norm(msh.inc_d[e]):
            raise SystemExit(f"隣接 {nodes[j]} の位置が member によって違う (group root {mem[0]})")
    groot = msh.root[np.array(nodes)]
    first = {}
    root_l = np.empty(n, dtype=np.int64)
    for i, r in enumerate(groot):
        root_l[i] = first.setdefault(int(r), i)
    root_l[:len(mem)] = 0
    sub = G.Mesh.__new__(G.Mesh)
    sub.nCells = n
    sub.inc_m = np.array([loc[int(x)] for x in msh.inc_m[es]], dtype=np.int64)
    sub.inc_j = np.array([loc[int(x)] for x in msh.inc_j[es]], dtype=np.int64)
    sub.inc_d0 = msh.inc_d[es].copy()
    sub.root = root_l
    sub.nmember = np.bincount(root_l, minlength=n)[root_l]
    sub._alpha = None
    return sub, P, nodes


def grad_at_root(sub, P, s):
    """形状固定で s 倍に縮小した stencil の LSQ 勾配 (lsq_merged_ref) と解析勾配 a の差 (成分ごと絶対値)。"""
    sub.inc_d = s * sub.inc_d0
    sub._alpha = None
    X = s * P
    phi = X @ A_LIN + 0.5 * np.einsum("ni,ij,nj->n", X, HESS, X)
    X = np.nan_to_num(X)
    g = G.lsq_merged_ref(sub, phi)                       # phi は double のまま渡す (縮小後の丸めを持ち込まない)
    return np.abs(g[0] - A_LIN), g[0]


def main():
    global SCR
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=SCR)
    ap.add_argument("--mesh", default="box32_j0.2")
    ap.add_argument("--extra-corner-meshes", default="box16_j0.2,box64_j0.2",
                    help="参考 (判定外): 他のジッタ格子の角 group")
    a = ap.parse_args()
    SCR = a.scratch
    md = os.path.join(SCR, "g_harness_mesh", a.mesh)
    msh = G.Mesh(os.path.join(md, "mesh.h5"), os.path.join(md, "bcondConfig.yaml"))
    h0 = 2 * np.pi / 32
    out = ["# G2′ 固定形状縮小試験 (plan boundary-node-periodic-gradient-fix §5.1 #8b、codex result-2 M2)",
           f"harness revision: {G.git_rev()}",
           f"メッシュ: {md}/mesh.h5 (32³、ジッタ ±0.2h、h = 2π/32 = {h0:.5f})",
           "参照: gharness.lsq_merged_ref (G2 の合併 stencil LSQ の double 参照、GPU と ≤ 1e-5·S で一致済み) をそのまま呼ぶ",
           f"場: φ = a·x + ½ xᵀH x、原点 = group の root、a = {tuple(A_LIN)}、H = {KAPPA:g}·[[1,0,−0.35],[0,0.5,0.15],[−0.35,0.15,0]]。"
           "解析勾配 (原点) = a。誤差 = 成分ごとの |∇_ref − a|",
           "縮小: 変位 d → s·d (s = 1, 1/2, 1/4、無次元形状固定)。h_s = s·h",
           "",
           "注記: 形状を固定した縮小では LSQ の二次場の誤差は s に厳密に比例する (M = Σ w d dᵀ は s に不変、右辺は s に比例) ので、"
           "次数はほぼ 1 になる。示すのは合併作用素の 1 次整合であり、ジッタ系列 (G2p_jitter.txt、64→128 の継ぎ目次数 0.855) の"
           "「判定不能」は診断としてそのまま残す (この試験で覆さない)。", ""]
    import collections
    cnt = collections.Counter(len(v) for v in msh.groups.values())
    out.append(f"実在の継ぎ目 group (member 数: 個数): {dict(sorted(cnt.items()))}")

    # 自己検査用: 全メッシュで同じ二次場 (実座標、root 原点ではなく大域座標) の lsq_merged_ref
    xg = msh.xyz.astype(np.float64)
    w, clean, _ = G.wrapped_coords(msh)
    phi_full = w @ A_LIN + 0.5 * np.einsum("ni,ij,nj->n", w, HESS, w)
    g_full = G.lsq_merged_ref(msh, phi_full)

    kinds = {2: "面 (2)", 4: "辺 (4)", 8: "角 (8)"}
    rows = []
    for nm, label in kinds.items():
        roots = sorted(r for r, v in msh.groups.items() if len(v) == nm and clean[r])
        pick = [roots[int(i)] for i in np.linspace(0, len(roots) - 1, min(N_PER_KIND, len(roots)))] if roots else []
        for r in dict.fromkeys(pick):
            rows.append((label, r, [r] + [m for m in msh.groups[r] if m != r]))
    interior = np.nonzero((msh.nmember == 1) & clean)[0]
    pick_i = [int(interior[int(i)]) for i in np.linspace(0, len(interior) - 1, N_INTERIOR)]
    for r in pick_i:
        rows.append(("内部 (1)", r, [r]))

    out.append("")
    out.append("## group ごと (誤差 e_x/e_y/e_z、次数は h→h/2 / h/2→h/4、C = e/h (h = h_s、最細))")
    out.append("| 種別 | root | member | 自己検査 max|Δ| (s=1、全メッシュとの差) | e (s=1) x / y / z | e (s=1/4) x / y / z | "
               "次数 h→h/2 x / y / z | 次数 h/2→h/4 x / y / z | 最小次数 | C = e/h (最細) x / y / z |")
    out.append("| --- " * 10 + "|")
    res = []
    for label, r, mem in rows:
        sub, P, nodes = sub_mesh(msh, mem)
        # 自己検査: s=1 で実座標の場 (w は wrapped 大域座標; 部分メッシュ上では root を原点へ平行移動した同じ一次・二次場)
        sub.inc_d = sub.inc_d0; sub._alpha = None
        wl = w[r] + P
        phil = wl @ A_LIN + 0.5 * np.einsum("ni,ij,nj->n", wl, HESS, wl)
        gs = G.lsq_merged_ref(sub, np.nan_to_num(phil))[0]
        selfchk = float(np.abs(gs - g_full[r]).max())
        E = []
        for s in SCALES:
            e, _ = grad_at_root(sub, P, s)
            E.append(e)
        E = np.array(E)
        o1 = np.log(E[0] / E[1]) / np.log(2.0)
        o2 = np.log(E[1] / E[2]) / np.log(2.0)
        C = E[2] / (h0 * SCALES[2])
        omin = float(np.nanmin(np.concatenate([o1, o2])))
        above = (E > DFLOOR).all(axis=0)                  # 全縮小率で丸めを超える成分
        omin_above = float(np.nanmin(np.concatenate([o1[above], o2[above]]))) if above.any() else float("nan")
        res.append(dict(label=label, r=r, nm=len(mem), omin=omin, omin_above=omin_above, below=[COMPS[c] for c in range(3) if not above[c]],
                        C=C, selfchk=selfchk, E=E))
        f3 = lambda v, fmt: " / ".join(format(x, fmt) for x in v)
        out.append(f"| {label} | {r} | {len(mem)} | {selfchk:.1e} | {f3(E[0], '.3e')} | {f3(E[2], '.3e')} | "
                   f"{f3(o1, '.4f')} | {f3(o2, '.4f')} | {omin:.4f} | {f3(C, '.4f')} |")

    seam = [x for x in res if x["nm"] > 1]
    inter = [x for x in res if x["nm"] == 1]
    out.append("")
    out.append("## 誤差定数 C = e/h (成分の最大) の比較")
    out.append("| 区分 | 個数 | C 最小 | C 中央値 | C 最大 |")
    out.append("| --- | --- | --- | --- | --- |")
    for lab in list(kinds.values()) + ["内部 (1)"]:
        cs = [x["C"].max() for x in res if x["label"] == lab]
        if cs:
            out.append(f"| {lab} | {len(cs)} | {min(cs):.4f} | {np.median(cs):.4f} | {max(cs):.4f} |")
    selfmax = max(x["selfchk"] for x in res)
    out.append(f"\n自己検査 (部分メッシュ抽出 = 全メッシュの lsq_merged_ref、s=1): 最大差 {selfmax:.1e} "
               f"({'一致' if selfmax <= 1e-12 else '不一致'}、判定 ≤ 1e-12)")

    # 参考 (判定外): 他のジッタ格子の角 group
    extra = []
    for em in [x for x in a.extra_corner_meshes.split(",") if x]:
        emd = os.path.join(SCR, "g_harness_mesh", em)
        if not os.path.exists(os.path.join(emd, "mesh.h5")):
            extra.append(f"- {em}: メッシュ無し")
            continue
        em_msh = G.Mesh(os.path.join(emd, "mesh.h5"), os.path.join(emd, "bcondConfig.yaml"))
        for r, v in em_msh.groups.items():
            if len(v) != 8:
                continue
            sub, P, _ = sub_mesh(em_msh, [r] + [m for m in v if m != r])
            E = np.array([grad_at_root(sub, P, s)[0] for s in SCALES])
            o = np.concatenate([np.log(E[0] / E[1]), np.log(E[1] / E[2])]) / np.log(2.0)
            hN = 2 * np.pi / round(2 * np.pi / np.median(np.linalg.norm(em_msh.inc_d, axis=1)))
            extra.append(f"- {em} 角 group root {r} (member 8): 最小次数 {np.nanmin(o):.4f}、e (s=1) = "
                         + " / ".join(f"{x:.3e}" for x in E[0]))
    if extra:
        out.append("\n## 参考 (判定外): 32³ 以外のジッタ格子の角 group (同じ生成器、別の格子)")
        out += extra

    # 判定
    n_by = {lab: sum(1 for x in seam if x["label"] == lab) for lab in kinds.values()}
    omin_seam = min(x["omin"] for x in seam)
    enough = all(n_by[lab] >= 2 for lab in kinds.values()) and len(seam) >= 6
    out.append("")
    out.append(f"継ぎ目 group の個数: " + "、".join(f"{k} {v}" for k, v in n_by.items()) + f"、計 {len(seam)}")
    out.append(f"継ぎ目の最小次数 (全成分・両対): {omin_seam:.4f}、内部の最小次数: {min(x['omin'] for x in inter):.4f}")
    belows = [f"{x['label']} root {x['r']} の {','.join(x['below'])} 成分 (e = " + " / ".join(f"{e:.1e}" for e in x['E'][:, COMPS.index(x['below'][0])]) + ")"
              for x in seam if x["below"]]
    if belows:
        out.append("誤差が全縮小率で倍精度丸め (≤ %.0e) 以下の成分 (次数を定義できない。§5.1 #8b の基準に床の扱いは無いので判定には使わず事実のみ): "
                   % DFLOOR + "; ".join(belows))
        out.append(f"それらを除いた継ぎ目の最小次数: {min(x['omin_above'] for x in seam):.4f}")
    if selfmax > 1e-12:
        verdict = "判定不能 (部分メッシュ抽出の自己検査が不一致)"
    elif not enough:
        short = [k for k in kinds.values() if n_by[k] < 2]
        verdict = (f"判定保留 (基準の個数条件を満たせない: {'、'.join(short)} の実在 group が 32³ ジッタ格子に "
                   f"{', '.join(str(cnt.get(nm, 0)) for nm, lab in kinds.items() if lab in short)} 個しか無い "
                   f"(三重周期立方体の 8 隅は 1 つの group)。取った group の最小次数 {omin_seam:.4f} は記録のみ)")
    else:
        verdict = ("PASS" if omin_seam >= 0.9 else "FAIL") + f" (継ぎ目の最小次数 {omin_seam:.4f} {'≥' if omin_seam >= 0.9 else '<'} 0.9)"
    out.append(f"\nVERDICT G2′ 固定形状 (面・辺・角 各 ≥ 2、計 ≥ 6、全成分の次数 ≥ 0.9): {verdict}")
    txt = "\n".join(out) + "\n"
    p = os.path.join(HERE, "G2p_fixedshape.txt")
    open(p, "w").write(txt)
    print(txt)
    print("->", p)


if __name__ == "__main__":
    main()
