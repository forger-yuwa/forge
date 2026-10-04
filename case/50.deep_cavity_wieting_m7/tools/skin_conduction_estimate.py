#!/usr/bin/env python3
"""薄肉過渡法の**横方向伝導**が、深部の見かけ熱流束をどれだけ作りうるか (W70 の計測モデル)。

W70 (TN D-5908) は後壁・床に 0.012 in (0.305 mm) の 304SS 薄板を貼り、
  q_meas(s) = ρ c τ (dT/dt)_s
で局所熱流束を出す。本文は「表面伝導・放射の補正はしていない」と明記している。
薄板の 1 次元エネルギー式は
  ρ c τ ∂T/∂t = q_conv(s) + λ_s τ ∂²T/∂s²      (背面は断熱シリカ)
なので、**上端が先に温まると、その後の (dT/dt) には横方向伝導で運ばれた分が混ざる**。
挿入 (0.8 s) の間に上端だけが ΔT 温まった状態で読むと、見かけの熱流束は

  q_apparent(s) ≈ q_conv(s) + λ_s τ ΔT / (L_hot · L_abs)

程度になる。ここで L_hot は上端の高熱流束帯の長さ、L_abs はそれを吸収する深さ幅。
"""
import numpy as np

LAM_S = 15.0        # 304SS の熱伝導率 [W/(m·K)] (300-400 K)
TAU = 0.305e-3      # 薄板厚 [m] (0.012 in)
RHO_C = 7900 * 500  # ρc [J/(m³·K)]
Q_FP = 29.3e3       # 文献の平板基準 [W/m²]
D = 20.32e-3        # キャビティ深さ
W = 1.270e-3        # 幅 (w/d = 0.063)

print("W70 の薄肉スキン横方向伝導の見積り")
print(f"  λ_s τ = {LAM_S*TAU:.4e} W/K   (1 m 幅あたりの伝導コンダクタンス)")
print()
print(f"{'ΔT_top [K]':>10} {'L_hot [mm]':>10} {'L_abs [mm]':>10} {'q_app [kW/m²]':>14} {'q_app/q_fp':>11}")
for dT in (5.0, 10.0, 22.0):                   # 本文: 計測区間の最大温度上昇 <22 K
    for L_hot in (1.0e-3, 2.0e-3):
        for L_abs in (2.0e-3, 5.0e-3):
            q_lat = LAM_S * TAU * dT / L_hot   # [W/m of width] 上端から流れ出す熱
            q_app = q_lat / L_abs              # それを L_abs の深さ帯で受けたときの見かけ面熱流束
            print(f"{dT:10.1f} {L_hot*1e3:10.1f} {L_abs*1e3:10.1f} {q_app*1e-3:14.2f} {q_app/Q_FP:11.4f}")
print()
print("対比 (NASA TN D-5908 Fig 6(a), w/d = 0.063 の実測):")
print("   x/d = 0.10 (= 1.6 すきま幅)  q_s/q_fp ≈ 0.17")
print("   x/d = 0.20                  q_s/q_fp ≈ 0.07")
print("   x/d ≥ 0.40                  q_s/q_fp ≈ 0.02 (計測限界帯)")
print()
print("→ 見積り値は実測の深部超過と**同じ桁**。W70 の深部値は横方向伝導の寄与を含みうる。")
print("   (本文自身が『表面伝導・放射の補正なし』と明記。Fig 10(a) でも実測は Burggraf 理論の上に乗る)")

# 特性長: 伝導が効く範囲
tau_ins = 0.8                                   # 挿入時間 [s]
alpha_s = LAM_S / RHO_C
L_diff = np.sqrt(alpha_s * tau_ins)
print(f"\n参考: 挿入 {tau_ins} s の間にスキン内を熱が拡散する長さ √(α_s t) = {L_diff*1e3:.2f} mm "
      f"(= {L_diff/W:.1f} すきま幅, x/d = {L_diff/D:.3f})")
