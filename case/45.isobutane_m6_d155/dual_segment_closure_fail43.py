"""plan axisymmetric-freestream-hoop-gauge §4.4 の記録: §4.3 の B の全域の判定で超えた CV の内訳。dual_segment_closure.py を同じプロセスで実行した後に exec する (名前空間を共有)。"""
# §4.3 の B の全域の判定で超えた CV の内訳 (記録用、dual_segment_closure.py の後に同じプロセスで exec する)
Sx = node_sum(exact_fx, bexact_x); Sy = node_sum(exact_f, bexact); Ax = node_abs(exact_fx, bexact_x); Ay = node_abs(exact_f, bexact)
ox = np.abs(Sx) > 100 * eps * (Ap + Ax); oy = np.abs(Sy - Ap) > 100 * eps * (Ap + Ay)
for nm, o, S_, A_, tgt in (("x", ox, Sx, Ax, 0.0), ("y", oy, Sy, Ay, Ap)):
    jj = j[o]
    print(f"{nm}: 超えた CV {o.sum()}、j の分布 {np.bincount(jj, minlength=121)[[0,1,2,3,4,5,60,110,115,117,118,119,120]].tolist()} (j 0,1,2,3,4,5,60,110,115,117,118,119,120)")
    ratio = np.abs(S_ - tgt)[o] / (eps * (Ap[o] + A_[o]))
    print(f"   |欠損|/(ε·(A+Σ|W|)) の中央値 {np.median(ratio):.0f}・最大 {ratio.max():.0f}")
    # 座標の桁落ちの目安: ε·|x|·(CV の周の長さの規模)/A  ではなく、ε·max|座標|/最短の辺 を相対誤差の目安とする
    nodes = np.flatnonzero(o)[:5]
    for n in nodes:
        print(f"   例 node {n}: j {j[n]} i {n // NJ} x/rt {x[n]/RT:.3f} r/rt {y[n]/RT:.4f} |欠損|/A {abs(S_[n]-(tgt[n] if np.ndim(tgt) else tgt))/Ap[n]:.2e} Σ|W|/A {A_[n]/Ap[n]:.2e}")
# 壁の第一層の厚さと座標の比
wall = np.flatnonzero(j == 120); fi = wall - 1
t = np.hypot(x[wall] - x[fi], y[wall] - y[fi])
print(f"壁と第一内部の節点の距離 / 座標の大きさ: 中央値 {np.median(t / np.hypot(x[wall], y[wall])):.2e}・最小 {np.min(t / np.hypot(x[wall], y[wall])):.2e}")
