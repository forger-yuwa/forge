帯修正 A/B (plan verification-m6-axis-wave-mesh-su2) の ParaView 用ファイル。座標は r_t 単位 (r_t = 0.07702 m)、軸対称の上半分。
  B1.xmf          方式 E の壁 (run_0046_ns_band_edge/res_12000)
  B0.xmf          adaptive の壁 (run_0048_ns_band_adaptive_ext/res_6000 = run_0045 から延長)
  A0p.xmf         旧最終壁 (run_0042_ns_restart_ctrl/res_12000、新バイナリ)
  Euler.xmf       設計壁の Euler (run_0047_euler_rt77p02_newbin)
  B1_minus_B0.xmf B1 − B0 (同じ (i, j) 節点)
変数: Mach, M_over_6_minus1_pct, M_vs_Euler_pct (Euler を (x, r/r_w) で補間して比), dpdx_nd = (r_t/p) ∂p/∂x, P, T, rho, Ux, Uy; 差: dM_pct, d_dpdx_nd
起動: ~/opt/ParaView-6.1.1-MPI-Linux-Python3.12-x86_64/bin/paraview --script=load_band_ab.py
生の結果 (全変数) は各 run の res_*.xmf をそのまま開いてもよい。
