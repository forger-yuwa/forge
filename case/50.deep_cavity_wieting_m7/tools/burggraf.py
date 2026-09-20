#!/usr/bin/env python3
"""Burggraf (1965) の非粘性コア理論による q_s/q_fp — NASA TN D-5908 Appendix B (B1)–(B4) の実装。

  q_s/q_fp = 0.6 Q(X) = 0.21/sqrt(1 + d/w) · [ ζ(1/2, s/(2(w+d))) − ζ(1/2, (s+w)/(2(w+d))) ]   (B4)
  ζ(1/2, v) = 1/sqrt(v) + Σ_{n=0..5} a_n (v+1)^n                                              (B3)
  適用条件: N*_Re,L > N*_Re,L,cr = 240 (L/w)^{4/3} (1 + d/w)                                    (B1)

s はキャビティ周長に沿う距離 (s=0 が後壁上端)。X = 1 − s/s0、s0 = 2d + w。
**理論は Pr=1 を仮定**しており、その下では壁温は q_s/q_fp の比から落ちる
(コアは全エンタルピー、平板の回復エンタルピーも Pr=1 なら全エンタルピー)。
実際は Pr<1 なので、比は (h_t − h_w)/(h_aw − h_w) の分だけ壁温に依存する。
"""
import numpy as np

A = [0.803323, -3.89728, 2.55002, -1.19121, 0.308284, -0.0335024]


def zeta_half(v):
    v = np.asarray(v, float)
    poly = sum(a * (v + 1.0) ** n for n, a in enumerate(A))
    return 1.0 / np.sqrt(v) + poly


def qs_over_qfp(s, w, d):
    """周長座標 s [m] における q_s/q_fp (B4)。"""
    s = np.asarray(s, float)
    k = 2.0 * (w + d)
    return 0.21 / np.sqrt(1.0 + d / w) * (zeta_half(s / k) - zeta_half((s + w) / k))


def re_crit(L, w, d):
    """(B1) 粘性コアになる臨界 Re* (Eckert 基準温度)。"""
    return 240.0 * (L / w) ** (4.0 / 3.0) * (1.0 + d / w)


if __name__ == "__main__":
    w, d = 1.270e-3, 20.32e-3
    L = 158.496e-3 - w
    print(f"w/d = {w/d:.3f}  (d/w = {d/w:.1f})")
    print(f"(B1) N*_Re,L,cr = {re_crit(L, w, d):.3e}   [L = {L*1e3:.1f} mm]")
    print(f"     参考: L = 88.9 mm (圧力孔位置) なら {re_crit(0.0889, w, d):.3e}")
    print(f"     W70 Table (本文): 2.50 × 10^n  → L = 157 mm 側が一致 (幾何の独立確認)")
    print(f"     試験 N*_Re,L = 1.4–1.8e5 < cr → **粘性コア** (理論は適用外, W70 本文と同じ)")
    print()
    print(f"{'x/d':>7} {'s [mm]':>8} {'理論 q_s/q_fp':>14}")
    for xd in (0.0125, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0):
        s = xd * d
        print(f"{xd:7.4f} {s*1e3:8.3f} {float(qs_over_qfp(s, w, d)):14.4f}")
