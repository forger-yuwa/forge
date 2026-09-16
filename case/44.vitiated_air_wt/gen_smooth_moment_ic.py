#!/usr/bin/env python3
"""凝縮モーメントの滑らかな時間次数試験用 IC (plan species-passive-scalar-unification §6-6, 2026-09-17)。
収束した凝縮流れの IC (run_0312 ic.h5: 節点座標 MESH/COORD, 保存量/原始量 VALUE) から液滴場を全部消し、湿り域 (S>1 の外側流線) の
中にガウス状の液滴ブロブ (g0 exp(−(x−xc)²/sx² − (y−yc)²/sy²), 液滴半径 r_d 一定, 実現可能域の内部 x=y=0.9) を置く。
核生成は solverConfig の condSigmaScale 3.0 で抑止する (J ∝ exp(−σ³…) で実効 0; 成長・蒸発は σ に依らない)。
使い方: gen_smooth_moment_ic.py SRC_ic.h5 DST_ic.h5 [--g0 1e-3] [--rd 1.5e-7] [--xy 0.9] [--mode wet|core]
  wet : 湿り域 (rog>1e-6) の x 中央 (55 %) と、その x での湿り域 y 中央に置く (run_0342–0345)
  core: 湿り域 x の 60 % 位置、壁半径の 45 % に置く (乾き核 → 蒸発消滅の試験; run_0338–0341)"""
import argparse, math, shutil
import h5py, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('src'); ap.add_argument('dst')
ap.add_argument('--g0', type=float, default=1.0e-3); ap.add_argument('--rd', type=float, default=1.5e-7); ap.add_argument('--xy', type=float, default=0.9)
ap.add_argument('--mode', choices=['wet', 'core'], default='wet')
ap.add_argument('--center', type=float, nargs=2, default=None, help='ブロブ中心 (x, y) を明示 (乾き IC など湿り域が無いとき)')
ap.add_argument('--sigma', type=float, nargs=2, default=None, help='ブロブ幅 (sx, sy)')
a = ap.parse_args()
with h5py.File(a.src, 'r') as h:
    X = np.asarray(h['MESH']['COORD']).reshape(-1, 3); V = h['VALUE']
    x, y = X[:, 0], X[:, 1]
    g = np.asarray(V['rog_0']) if 'rog_0' in V else np.zeros_like(x); T = np.asarray(V['T']); ro = np.asarray(V['ro'])
wet = g > 1e-6
if a.center is not None:
    xc, yc = a.center; sx, sy = a.sigma if a.sigma else (0.2, 0.07)
elif a.mode == 'wet':
    xc = x[wet].min() + 0.55*(x[wet].max() - x[wet].min()); sel = wet & (np.abs(x - xc) < 0.1)
    yc = 0.5*(y[sel].min() + y[sel].max()); sy = 0.35*(y[sel].max() - y[sel].min()); sx = 0.12*(x[wet].max() - x[wet].min())
else:
    xc = x[wet].min() + 0.6*(x[wet].max() - x[wet].min()); ywall = y[np.abs(x - xc) < 0.02].max()
    yc = 0.45*ywall; sy = 0.25*ywall; sx = 0.15*(x[wet].max() - x[wet].min())
rl = np.maximum(1000.0 - 0.12*(277.0 - T), 920.0)   # solver と同じ H2O 液密度
gb = a.g0*np.exp(-((x - xc)**2/sx**2 + (y - yc)**2/sy**2)); gb[gb < 1e-12] = 0.0
Q0 = gb/((4.0/3.0)*math.pi*rl*a.rd**3); Q1 = a.xy*Q0*a.rd; Q2 = a.xy*Q0*a.rd*a.rd
shutil.copy(a.src, a.dst)
with h5py.File(a.dst, 'r+') as o:
    W = o['VALUE']
    for k, v in (('g_0', gb), ('Q0_0', Q0), ('Q1_0', Q1), ('Q2_0', Q2), ('rog_0', ro*gb), ('roQ0_0', ro*Q0), ('roQ1_0', ro*Q1), ('roQ2_0', ro*Q2)):
        if k in W: W[k][...] = v.astype(W[k].dtype)
        else: W.create_dataset(k, data=v.astype(W['ro'].dtype))   # 乾き IC (モーメント無し) には新設
    for k in ('condClampCorrQ_0', 'condClampCorr_0', 'condLim_0', 'condR30_0', 'condDrdt_0', 'condTheta_0'):
        if k in W: W[k][...] = 0*np.asarray(W[k])
print(f'blob center ({xc:.4f}, {yc:.4f}) sigma ({sx:.4f}, {sy:.4f}) r_d {a.rd:g} x=y={a.xy}: cells {int((gb > 0).sum())} (in wet zone {int((wet & (gb > 1e-6)).sum())} of {int((gb > 1e-6).sum())}), max g {gb.max():.3e}, max Q0 {Q0.max():.3e}')
