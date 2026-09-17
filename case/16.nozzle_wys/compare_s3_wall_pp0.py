#!/usr/bin/env python3
"""S3 (speciesFaceReconstruction 2) 有無で Wyslouzil Fig.3 の壁圧 p/p0 を比較する
   (plan species-passive-scalar-unification §6-5; run_0482 [S3] vs run_0483 [旧経路] vs 実験)。
usage: compare_s3_wall_pp0.py [OUT.png]"""
import sys, glob
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
for fp in glob.glob('/home/sano/.fonts/*Noto*CJK*'):
    font_manager.fontManager.addfont(fp)
plt.rcParams['font.family'] = 'Noto Sans CJK JP'
plt.rcParams['axes.unicode_minus'] = False

S = 'wall_pp0_csv/'
exp = np.array([[float(v) for v in l.split(',')[:3]] for l in open('wyslouzil_fig3_pp0.csv').read().strip().split('\n')[1:]])
xe, iso, cond = exp[:, 0]*10, exp[:, 1], exp[:, 2]        # cm -> mm
a = np.genfromtxt(S + 'wall_run_0482_passive_wys_s1_sfr2_c1.csv', delimiter=',', names=True)   # S3
b = np.genfromtxt(S + 'wall_run_0483_passive_wys_s0_sfr0.csv',   delimiter=',', names=True)   # 旧経路
out = sys.argv[1] if len(sys.argv) > 1 else 's3_vs_legacy_wall_pp0.png'

def onset(d, thr=1e-3):
    m = d['center_g'] > thr
    return d['x_mm'][m][0] if m.any() else float('nan')

fig = plt.figure(figsize=(13.2, 5.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 0.72], wspace=0.20, hspace=0.42)
ax = fig.add_subplot(gs[:, 0]); az = fig.add_subplot(gs[0, 1]); ad = fig.add_subplot(gs[1, 1])
oa, ob = onset(a), onset(b)
for A in (ax, az):
    A.plot(xe, iso, 's', ms=5.5, mfc='none', mew=1.1, color='#7f8c8d', label='実験 等エントロピー (dry)')
    A.plot(xe, cond, 'o', ms=7, mfc='none', mew=1.8, color='k', label='実験 Wyslouzil Fig.3 (1.00 kPa, 凝縮)')
    A.plot(b['x_mm'], b['contour_wall'], color='#1f4e79', lw=2.1, label='forge 旧経路 (SFR 0)')
    A.plot(a['x_mm'], a['contour_wall'], color='#c0392b', lw=2.1, ls='--', label='forge S3 (SFR 2 + coupling 1)')
    A.grid(alpha=.3)
ax.set_xlabel('軸方向位置 x [mm]  (スロート = 0)'); ax.set_ylabel('壁静圧 p / p0     (p0 = 59.07 kPa)')
ax.set_xlim(-10, 95); ax.set_ylim(0.15, 0.78)
ax.set_title('Wyslouzil ノズルの壁圧分布', loc='left')
ax.axvspan(11.2, 21.4, color='#f0c419', alpha=.16)
ax.annotate('実験はこの区間で等エントロピーから外れる\n(測点間隔 10 mm なので onset は 11–21 mm)',
            xy=(16.3, 0.335), xytext=(24, 0.56), fontsize=9, color='#8a6d00',
            arrowprops=dict(arrowstyle='->', color='#8a6d00', lw=1.1))
ax.legend(fontsize=9, framealpha=.92, loc='upper right')

az.set_xlim(14, 46); az.set_ylim(0.325, 0.372)
az.set_title('凝縮バンプの拡大', loc='left', fontsize=10)
for x0, c, lab, dy in ((ob, '#1f4e79', f'旧経路 onset {ob:.2f}', 0.3695), (oa, '#c0392b', f'S3 onset {oa:.2f}', 0.3665)):
    az.axvline(x0, color=c, lw=1.1, ls=':')
    az.annotate(lab + ' mm', xy=(x0 + 0.4, dy), color=c, fontsize=8.8, va='center')

ad.plot(a['x_mm'], (a['contour_wall'] - b['contour_wall'])*1e3, color='#c0392b', lw=1.8)
ad.axhline(0, color='k', lw=.8); ad.grid(alpha=.3)
ad.set_xlim(-10, 95); ad.set_xlabel('軸方向位置 x [mm]')
ad.set_ylabel('S3 − 旧経路\n[p/p0, x1e-3]', fontsize=9)
d = a['contour_wall'] - b['contour_wall']; k = int(np.argmax(np.abs(d)))
ad.annotate(f'最大 {d[k]*1e3:+.2f}e-3  ({abs(d[k])/b["contour_wall"][k]*100:.2f} %)  @ x {a["x_mm"][k]:.0f} mm',
            xy=(a['x_mm'][k], d[k]*1e3), xytext=(a['x_mm'][k] + 8, d[k]*1e3*0.75),
            fontsize=8.8, color='#c0392b', arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.0))
fig.suptitle('スカラ移流を 2 次 (S3) にすると凝縮 onset が +0.85 mm 下流へ動く — case/16 run_0483 (SFR 0) vs run_0482 (SFR 2)', fontsize=12)
fig.subplots_adjust(top=0.90, left=0.075, right=0.985, bottom=0.10)
fig.savefig(out, dpi=140)
print('wrote', out)
print(f'onset (中心線 g>1e-3): 旧経路 {ob:.2f} mm, S3 {oa:.2f} mm, 差 {oa-ob:+.2f} mm')
d = np.abs(a['contour_wall'] - b['contour_wall']); i = int(np.argmax(d))
print(f'壁圧の最大差 {d[i]:.5f} (相対 {d[i]/b["contour_wall"][i]*100:.2f} %) at x = {a["x_mm"][i]:.2f} mm')
for lab, dd in (('旧経路', b), ('S3', a)):
    r = np.interp(xe, dd['x_mm'], dd['contour_wall'])
    print(f'  実験 (凝縮) との差: {lab} 平均 {np.mean(r-cond):+.4f}, 最大 |差| {np.max(np.abs(r-cond)):.4f}')
