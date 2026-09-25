#!/usr/bin/env python3
"""G0 の CPU 作用素試験 (plan boundary-node-periodic-gradient-fix §6 G0)。GPU でメッシュを作らない 3 ケース。

`calcGradient_d.cu` lsqPre_mergePeriodic の式を CPU に移した `merged_coefs` で係数を焼き (double → float32)、
実行時と同じく float32 で Σ c·Δφ を積む。

1. float32 反例 (codex 実装レビュー M1): 並進 (100,0,0) の継ぎ目で、root 側 incidence の変位
   d0 = f32(0.03) − f32(0)、member 側 d1 = f32(100.03) − f32(100)。同じ物理隣接 (同 root, |d0−d1| ≤ 1e-4 h_min)
   なので同値類 (α=1/2)。plan の合格: 1.000000 (root 両順序)。
   - 模型 A (plan §4.1 の想定): 各部分 CV が読む Δφ_i = φ(j_i) − φ(m_i) が自分の格納座標の線形場 (Δφ_i = d_i)。
   - 模型 B (参考): DOF 共有の実態。j0≡j1、m0≡m1 なので Δφ は 1 つ (= 真の変位 0.03 を f32 で格納した値の差)。
     座標の丸め (|x|=100 で ulp 7.6e-6、h=0.03 の 2.5e-4) が入力側の床として残る。
   比較として「代表変位 (同値類の最初の incidence の d) で組む旧実装」も出す。
2. 部分 CV で rank 欠損 → 合併で回復: root 側の部分は xy 平面内の隣接だけ (rank 2)、member 側は yz 平面内だけ (rank 2)。
   合併 stencil は rank 3。線形場で厳密に戻ること、部分ごとに解いて和を取る旧方式 (打ち切り込み) が戻らないことを示す。
3. 合併しても退化する (全隣接が xy 平面内) → 退化方向は打ち切りの参照解 (double、同じ打ち切り) と比較。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G   # noqa: E402

f32 = np.float32


def merged_coefs(X, root_first_members, incs, rootof, thresh=G.LSQ_THRESH, representative=False):
    """lsqPre_mergePeriodic の 1 group ぶん。X: 格納座標 (float32, (n,3))。root_first_members: [root, members...]。
    incs: [(m, j), ...] (m は group の partial、CSR 順)。rootof: 節点 → root。戻り値: 各 incidence の float32 係数 (k,3)。"""
    order = [e for m in root_first_members for e in incs if e[0] == m]
    d = np.array([X[j].astype(np.float64) - X[m].astype(np.float64) for m, j in order])
    L = np.linalg.norm(d, axis=1)
    tol = 1e-4 * L[L > 0].min()
    cls, lab = [], []
    for k, (m, j) in enumerate(order):
        found = -1
        for c, (jr, dd, cnt) in enumerate(cls):
            if jr == rootof[j] and np.linalg.norm(d[k] - dd) <= tol:
                found = c; break
        if found < 0:
            cls.append([rootof[j], d[k], 0]); found = len(cls) - 1
        cls[found][2] += 1; lab.append(found)
    dd = np.array([cls[lab[k]][1] for k in range(len(order))]) if representative else d
    cnt = np.array([cls[lab[k]][2] for k in range(len(order))], dtype=np.float64)
    w = 1.0 / np.sum(dd * dd, axis=1) / cnt
    M = np.einsum("k,ki,kj->ij", w, dd, dd)
    Minv, degen = G.pinv_sym3(M[None])
    c = (Minv[0] @ (w[:, None] * dd).T).T.astype(f32)
    return order, c, bool(degen[0]), M


def apply(order, c, dphi):
    g = np.zeros(3, dtype=f32)
    for k in range(len(order)):
        g = (g + c[k] * f32(dphi[k])).astype(f32)
    return g.astype(np.float64)


def case_counterexample(out):
    out.append("## 1. float32 反例 (d0 = f32(0.03)−f32(0), d1 = f32(100.03)−f32(100))")
    # 節点: 0=m (x=0), 1=m' (x=100), 2=j (x=0.03), 3=j' (x=100.03), 4,5 = root 側の y・z 隣接 (rank 補完)
    X = np.array([[0, 0, 0], [100, 0, 0], [0.03, 0, 0], [100.03, 0, 0], [0, 0.05, 0], [0, 0, 0.05]], dtype=f32)
    d0 = f32(0.03) - f32(0.0); d1 = f32(100.03) - f32(100.0)
    out.append(f"d0 = {float(d0):.10g}, d1 = {float(d1):.10g}, 相対差 {abs(float(d1) / float(d0) - 1):.2e} (照合許容 1e-4 h_min 内: 同値類)")
    out.append("| 模型 | root | 実装 | ∂φ/∂x | |1−g| | 判定 (1.000000) |")
    out.append("| --- | --- | --- | --- | --- | --- |")
    ok_all = True
    for rootfirst in (0, 1):
        m_r, m_m = (0, 1) if rootfirst == 0 else (1, 0)
        rootof = {0: min(0, 1), 1: min(0, 1), 2: 2, 3: 2, 4: 4, 5: 5}
        # root を入れ替える = group の index 最小がどちら側か。rootof の値は識別子なので同じでよい
        incs = [(0, 2), (0, 4), (0, 5), (1, 3)]
        for rep in (False, True):
            order, c, _, _ = merged_coefs(X, [m_r, m_m], incs, rootof, representative=rep)
            for model in ("A", "B"):
                if model == "A":        # 各部分が自分の格納座標の線形場を読む
                    dphi = [float(X[j][0]) - float(X[m][0]) for m, j in order]
                else:                   # DOF 共有: Δφ は 1 つ (真の変位 0.03 の f32 値)
                    dphi = [float(f32(0.03)) if X[j][0] - X[m][0] > 0.02 else float(X[j][0] - X[m][0])
                            for m, j in order]
                g = apply(order, c, dphi)
                e = abs(1.0 - g[0])
                txt = f"{g[0]:.6f}"
                judge = ""
                if model == "A" and not rep:
                    ok = txt == "1.000000"
                    ok_all &= ok
                    judge = "ok" if ok else "NG"
                out.append(f"| {model} | {'x=0 側' if rootfirst == 0 else 'x=100 側'} | {'代表変位 (旧)' if rep else '実変位 (現行)'} "
                           f"| {txt} | {e:.1e} | {judge or '参考'} |")
    out.append(f"VERDICT G0 float32 反例: {'PASS' if ok_all else 'FAIL'} (判定は模型 A・実変位の 2 行)")
    return ok_all


def case_rank_recovery(out):
    out.append("\n## 2. 部分 CV の rank 欠損 → 合併で回復")
    h = 0.1
    # root partial 0: 隣接は xy 平面内だけ (rank 2)。member partial 1 (並進 (10,0,0) の像): 隣接は yz 平面内だけ (rank 2)。
    # 両部分の張る空間は y で重なる → 部分ごとに解いて和を取る旧方式は y 成分を 2 重に数える
    X = np.array([[0, 0, 0], [10, 0, 0],
                  [h, 0, 0], [-0.5 * h, 0.8 * h, 0],                  # root 側 (非対称)
                  [10, -h, 0], [10, 0.3 * h, 1.1 * h], [10, 0.2 * h, -0.9 * h]], dtype=f32)
    rootof = {0: 0, 1: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    incs = [(0, 2), (0, 3), (1, 4), (1, 5), (1, 6)]
    a = np.array([0.7, -1.3, 0.4])
    # 線形場は各部分の局所座標 (member は並進を戻す) で与える: 実 DOF でも Δφ = a·d
    order, c, degen, M = merged_coefs(X, [0, 1], incs, rootof)
    dphi = [a @ (X[j].astype(np.float64) - X[m].astype(np.float64)) for m, j in order]
    g = apply(order, c, dphi)
    e = np.abs(g - a).max() / np.linalg.norm(a)
    ev = np.linalg.eigvalsh(M)
    out.append(f"合併 M の固有値比 λmin/λmax = {ev[0] / ev[2]:.3f} (打ち切り {degen})。線形場 a = {tuple(a)}")
    out.append(f"合併 stencil: g = ({g[0]:.7f}, {g[1]:.7f}, {g[2]:.7f})、最大誤差/|a| = {e:.2e}")
    # 旧方式: 部分ごとに解いて和 (各部分は自分の LSQ + 打ち切り)
    gs = np.zeros(3)
    for part in (0, 1):
        oi = [(m, j) for m, j in incs if m == part]
        o2, c2, dg, _ = merged_coefs(X, [part], oi, {k: k for k in range(7)})
        gs += apply(o2, c2, [a @ (X[j].astype(np.float64) - X[m].astype(np.float64)) for m, j in o2])
    e_old = np.abs(gs - a).max() / np.linalg.norm(a)
    out.append(f"旧方式 (部分ごとの LSQ + 打ち切り、の和): g = ({gs[0]:.4f}, {gs[1]:.4f}, {gs[2]:.4f})、最大誤差/|a| = {e_old:.2e} "
               "(重なる y 方向を 2 重に数える)")
    ok = e <= 1e-5 and e_old > 1e-2
    out.append(f"VERDICT G0 rank 回復: {'PASS' if ok else 'FAIL'} (閾値 1e-5)")
    return ok


def case_degenerate(out):
    out.append("\n## 3. 合併しても退化 (全隣接が xy 平面内) → 打ち切りの参照解と比較")
    h = 0.1
    X = np.array([[0, 0, 0], [10, 0, 0], [h, 0.2 * h, 0], [-0.8 * h, 0.1 * h, 0],
                  [10, 1.1 * h, 0], [10, -0.9 * h, 0], [10.05, 0.4 * h, 0]], dtype=f32)
    rootof = {k: k for k in range(7)}; rootof[1] = 0
    incs = [(0, 2), (0, 3), (1, 4), (1, 5), (1, 6)]
    a = np.array([0.7, -1.3, 0.4])            # z 成分は stencil から決まらない (退化方向)
    order, c, degen, M = merged_coefs(X, [0, 1], incs, rootof)
    d = np.array([X[j].astype(np.float64) - X[m].astype(np.float64) for m, j in order])
    dphi = d @ a
    g = apply(order, c, dphi)
    # 参照: double の M⁺ (同じ打ち切り) で解く
    w = 1.0 / np.sum(d * d, axis=1)
    Minv, _ = G.pinv_sym3((np.einsum("k,ki,kj->ij", w, d, d))[None])
    gref = Minv[0] @ (w[:, None] * d).T @ dphi
    e = np.abs(g - gref).max() / np.linalg.norm(gref)
    out.append(f"打ち切り {degen}。g = ({g[0]:.7f}, {g[1]:.7f}, {g[2]:.3e})、参照 = ({gref[0]:.7f}, {gref[1]:.7f}, {gref[2]:.3e})")
    out.append(f"非退化 (xy) 成分は a と一致: 最大誤差/|a_xy| = {np.abs(g[:2] - a[:2]).max() / np.linalg.norm(a[:2]):.2e}、"
               f"退化 (z) 成分 = {g[2]:.2e} (打ち切りで 0)")
    ok = e <= 1e-5 and np.abs(g[:2] - a[:2]).max() / np.linalg.norm(a[:2]) <= 1e-5
    out.append(f"VERDICT G0 退化方向: {'PASS' if ok else 'FAIL'} (参照解との差/|g_ref| = {e:.2e} ≤ 1e-5)")
    return ok


def main():
    out = ["# G0 CPU 作用素試験 (plan boundary-node-periodic-gradient-fix §6 G0)",
           f"harness revision: {G.git_rev()}  (lsqPre_mergePeriodic = calcGradient_d.cu:731-805 の式を CPU に移植)", ""]
    ok = case_counterexample(out)
    ok &= case_rank_recovery(out)
    ok &= case_degenerate(out)
    txt = "\n".join(out) + "\n"
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "G0_cpu_cases.txt")
    open(p, "w").write(txt)
    print(txt)
    print("->", p)


if __name__ == "__main__":
    main()
