#!/usr/bin/env python3
"""⑤ SERN 2 バンド構造メッシュ (スリットカウル) の整合テスト (S2)。"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec  # noqa: E402
from forge_design.meshing.mesh_sern import PHYS_SERN, SernMeshParams, generate_sern_mesh, write_msh41_named  # noqa: E402

FAIL = 0
def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1

k = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0), L_cowl=1.0,
                             x_max=6.0, nj=201, dx=4e-3, p_ext_over_p_in=0.05)).march()
d = k.design_ramp(M_c=3.3, f=0.45)
prm = SernMeshParams(ni_up=8, ni_noz=40, ni_plume=60, nj_top=31, nj_bot=21, scale=0.1, interface_angle=float(k.TH[-1, 0]),
                     top_ext_angle=d.info["theta_e"])
coords, quads, bedges, info, y_mid, y_top = generate_sern_mesh(d, prm)
ni, njt, njb, ite = info["ni"], info["nj_top"], info["nj_bot"], info["i_te"]
check("セル数 = (ni−1)(nj_top−1 + nj_bot−1)", quads.shape[0] == (ni - 1) * (njt - 1 + njb - 1))
check("ノード数 = ni·nj_bot + ni·(nj_top−1) + 重複 station 数", coords.shape[0] == ni * njb + ni * (njt - 1) + len(info["dup_stations"]))
# 反時計回り (符号付き面積 > 0)
xy = coords[:, :2]
a = np.zeros(len(quads))
for kk in range(4):
    p0, p1 = xy[quads[:, kk]], xy[quads[:, (kk + 1) % 4]]
    a += p0[:, 0] * p1[:, 1] - p1[:, 0] * p0[:, 1]
check("全 quad が反時計回り・非退化", np.all(a > 1e-14), f"min area {a.min():.2e}")
# 境界辺の数
check("cowl_in と cowl_out の辺数 = i_te", len(bedges["cowl_in"]) == ite and len(bedges["cowl_out"]) == ite)
check("ramp + top_out の辺数 = ni−1", len(bedges["ramp"]) + len(bedges["top_out"]) == ni - 1)
check("outlet 辺数 = nj_top−1 + nj_bot−1", len(bedges["outlet"]) == njt - 1 + njb - 1)
# スリット: cowl_in/cowl_out のノードは座標一致・ID 相異
ci = np.array(bedges["cowl_in"]); co = np.array(bedges["cowl_out"])
check("スリット: 上下ノードの座標一致", np.allclose(coords[ci], coords[co]))
check("スリット: 上下ノード ID は別 (TE 点のみ共有)", not np.any(ci[:, 0] == co[:, 0]) and ci[-1, 1] == co[-1, 1])
# 全境界辺の総和 = 閉境界 (各セル辺で共有されない辺の数)
from collections import Counter
cnt = Counter()
for q in quads:
    for kk in range(4):
        e = tuple(sorted((int(q[kk]), int(q[(kk + 1) % 4]))))
        cnt[e] += 1
free = {e for e, c in cnt.items() if c == 1}
tagged = {tuple(sorted(map(int, e))) for nm in bedges for e in bedges[nm]}
check("境界辺タグの集合 = quad の非共有辺の集合", free == tagged, f"free {len(free)} tagged {len(tagged)}")
# ランプ輪郭がメッシュ上境界に一致
rx = coords[[e[0] for e in bedges["ramp"]], 0] / prm.scale
ry = coords[[e[0] for e in bedges["ramp"]], 1] / prm.scale
m = rx >= 0
check("ランプ上境界 = 設計輪郭 (補間 1e-9)", np.allclose(ry[m], np.interp(rx[m], d.ramp_xy[:, 0], d.ramp_xy[:, 1]), atol=1e-9))
out = Path(__file__).with_name("_sern_mesh_test.msh"); write_msh41_named(out, coords, quads, bedges, PHYS_SERN)
txt = out.read_text(); check("msh4.1 書き出し (PhysicalNames 9 個)", "$PhysicalNames\n9" in txt); out.unlink()

# ---------- ext_top (ランプ側外部流ブロック, plan §4.11) ----------
prm2 = SernMeshParams(**{**prm.__dict__, "ext_top": True, "nj_ext_top": 13, "nj_wake": 5, "top_depth": 1.5})
c2, q2, b2, i2, ym2, yt2 = generate_sern_mesh(d, prm2)
ni2, njT, njW, kte = i2["ni"], i2["nj_ext_top"], i2["nj_wake"], i2["i_ramp_te"]
check("ext_top: wake ブロックあり (θ_e<0 なので base は有限)", i2["has_wake"] and i2["h_base"] > 0, f"h_base {i2['h_base']:.4f}")
n_old = (ni2 - 1) * (njt - 1 + njb - 1)
check("ext_top: セル数 = 旧 + top + wake", q2.shape[0] == n_old + (ni2 - 1) * (njT - 1) + (ni2 - 1 - kte) * (njW - 1))
check("ext_top: 先頭のセル・ノードは旧メッシュと同一 (順序互換)", np.array_equal(q2[:n_old], quads) and np.allclose(c2[:coords.shape[0]], coords))
a2 = np.zeros(len(q2)); xy2 = c2[:, :2]
for kk in range(4):
    p0, p1 = xy2[q2[:, kk]], xy2[q2[:, (kk + 1) % 4]]
    a2 += p0[:, 0] * p1[:, 1] - p1[:, 0] * p0[:, 1]
check("ext_top: 全 quad が反時計回り・非退化", np.all(a2 > 1e-14), f"min area {a2.min():.2e}")
cnt2 = Counter()
for q in q2:
    for kk in range(4):
        cnt2[tuple(sorted((int(q[kk]), int(q[(kk + 1) % 4]))))] += 1
free2 = {e for e, c in cnt2.items() if c == 1}; tagged2 = {tuple(sorted(map(int, e))) for nm in b2 for e in b2[nm]}
check("ext_top: 境界辺タグの集合 = 非共有辺の集合", free2 == tagged2, f"free {len(free2)} tagged {len(tagged2)}")
check("ext_top: vehicle 辺数 = 上面 k + base (nj_wake−1)", len(b2["vehicle"]) == kte + (njW - 1))
check("ext_top: ramp 辺数 = k (後縁まで), top_out 辺数 = ni−1", len(b2["ramp"]) == kte and len(b2["top_out"]) == ni2 - 1)
check("ext_top: inlet_ext 辺数 = nj_bot−1 + nj_ext_top−1", len(b2["inlet_ext"]) == njb - 1 + njT - 1)
check("ext_top: outlet 辺数 = 4 ブロック分", len(b2["outlet"]) == (njt - 1) + (njb - 1) + (njT - 1) + (njW - 1))
vy = c2[[e[0] for e in b2["vehicle"][:kte]], 1] / prm.scale
check("ext_top: 機体上面はランプより上", np.all(vy > d.ramp_xy[:, 1].max() - 1e-12))
check("ext_top: 燃焼器 IC 領域 (y_mid < y < y_top) はランプ上で閉じる", float(yt2(0.5 * d.L_ramp)) < i2["y_veh"])
out = Path(__file__).with_name("_sern_mesh_test_ext.msh"); write_msh41_named(out, c2, q2, b2, PHYS_SERN)
check("ext_top: msh4.1 書き出し (PhysicalNames 10 個)", "$PhysicalNames\n10" in out.read_text()); out.unlink()

# ---------- ext_top + vehicle_taper (base 無し、TE 共有) ----------
prm3 = SernMeshParams(**{**prm2.__dict__, "vehicle_taper": 1.0})
c3, q3, b3, i3, ym3, yt3 = generate_sern_mesh(d, prm3)
kte3 = i3["i_ramp_te"]
check("taper: wake 無し (h_base = 0)", (not i3["has_wake"]) and i3["h_base"] == 0.0)
check("taper: vehicle 辺数 = k (上面のみ)", len(b3["vehicle"]) == kte3)
te_ramp = b3["ramp"][-1][1]; te_veh = b3["vehicle"][-1][1]
check("taper: 後縁ノードをランプと機体上面で共有", te_ramp == te_veh)
vx = c3[[e[0] for e in b3["vehicle"]], 0] / prm.scale; vy = c3[[e[0] for e in b3["vehicle"]], 1] / prm.scale
ry3 = np.interp(vx, d.ramp_xy[:, 0], d.ramp_xy[:, 1]); ry3[vx < 0] = 1.0
check("taper: 機体上面はランプ以上 (厚さ ≥ 0)", np.all(vy >= ry3 - 1e-12))
a3 = np.zeros(len(q3)); xy3 = c3[:, :2]
for kk in range(4):
    p0, p1 = xy3[q3[:, kk]], xy3[q3[:, (kk + 1) % 4]]
    a3 += p0[:, 0] * p1[:, 1] - p1[:, 0] * p0[:, 1]
check("taper: 全 quad が反時計回り・非退化", np.all(a3 > 1e-14), f"min area {a3.min():.2e}")
cnt3 = Counter()
for q in q3:
    for kk in range(4):
        cnt3[tuple(sorted((int(q[kk]), int(q[(kk + 1) % 4]))))] += 1
free3 = {e for e, c in cnt3.items() if c == 1}; tagged3 = {tuple(sorted(map(int, e))) for nm in b3 for e in b3[nm]}
check("taper: 境界辺タグの集合 = 非共有辺の集合", free3 == tagged3, f"free {len(free3)} tagged {len(tagged3)}")
print("ALL PASS" if FAIL == 0 else f"{FAIL} FAIL"); sys.exit(1 if FAIL else 0)
