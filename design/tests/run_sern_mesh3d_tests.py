#!/usr/bin/env python3
"""⑤ SERN 3D hex メッシュ (mesh_sern3d) の整合テスト: 境界の閉性・タグ分割 (R2 の vehicle タグ)・重複ノード。"""
import sys
from collections import Counter
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec  # noqa: E402
from forge_design.meshing.mesh_sern3d import PHYS_SERN3D, SernMesh3DParams, generate_sern_mesh3d  # noqa: E402

FAIL = 0
def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1

k = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0), L_cowl=1.0,
                             x_max=6.0, nj=201, dx=4e-3, p_ext_over_p_in=0.05)).march()
d = k.design_ramp(M_c=3.3, f=0.45)
HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def closure(hexes, B):
    """各 hex 面: 2 hex で共有 or ちょうど 1 つの境界群に属する。"""
    cnt = Counter()
    for h in hexes:
        for f in HEX_FACES:
            cnt[frozenset(h[i] for i in f)] += 1
    bnd = Counter()
    for nm, faces in B.items():
        for q in faces:
            bnd[frozenset(q)] += 1
    unshared = [f for f, c in cnt.items() if c == 1]
    missing = sum(1 for f in unshared if bnd.get(f, 0) != 1)
    extra = sum(1 for f in bnd if cnt.get(f, 0) != 1)
    return len(unshared), missing, extra


for label, prm in (("W_vehicle=None (旧: 遠方境界まで機体)", SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6,
                                                                 W=2.0, Z_ext=1.5, interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"])),
                   ("W_vehicle=3.0", SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6,
                                                      W=2.0, Z_ext=1.5, W_vehicle=3.0, interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"])),
                   ("外側空間なし (nz_out=0)", SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=0,
                                                         W=2.0, interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"]))):
    coords, hexes, B, info, y_mid = generate_sern_mesh3d(d, prm)
    print(f"--- {label}: cells {info['cells']} nodes {info['nodes']} vehicle faces {info['n_vehicle_faces']}")
    nun, miss, extra = closure(hexes, B)
    check("境界の閉性: 非共有 hex 面 = 境界 quad (漏れ 0, 余り 0)", miss == 0 and extra == 0, f"unshared {nun} missing {miss} extra {extra}")
    top_all = len(B["ramp"]) + len(B["vehicle"]) + len(B["underside_far"]) + len(B["top_out"])
    check("上面 quad 数 = (ni−1)(nz−1)", top_all == (info["ni"] - 1) * (info["nz"] - 1))
    i_ramp = int(np.sum(0.5 * (np.r_[coords[:0, 0]] if False else 0)))  # placeholder (未使用)
    # ramp はノズル幅内 (z ≤ W/2) だけ
    zc = lambda faces: np.array([coords[list(q), 2].mean() for q in faces]) / prm.scale
    if len(B["ramp"]):
        check("ramp の面重心 z ≤ W/2", np.all(zc(B["ramp"]) <= 0.5 * prm.W + 1e-12), f"max z {zc(B['ramp']).max():.3f}")
    if prm.nz_out > 0:
        check("vehicle の面重心 z > W/2", np.all(zc(B["vehicle"]) > 0.5 * prm.W - 1e-12))
        if prm.W_vehicle is not None:
            check("vehicle の面重心 z ≤ W_vehicle/2、超えた分は underside_far", np.all(zc(B["vehicle"]) <= 0.5 * prm.W_vehicle + 1e-12)
                  and len(B["underside_far"]) > 0 and np.all(zc(B["underside_far"]) > 0.5 * prm.W_vehicle - 1e-12), f"max vehicle z {zc(B['vehicle']).max():.3f}")
        else:
            xt = np.array([coords[list(q), 0].mean() for q in B["top_out"]]) / prm.scale
            check("W_vehicle=None: x ≤ L_ramp の上面は ramp/vehicle だけ (top_out は後縁より下流のみ)", np.all(xt > info["L_ramp"]))
    else:
        check("外側空間なし: vehicle 面は 0", len(B["vehicle"]) == 0)
    check("physID は一意", len(set(PHYS_SERN3D.values())) == len(PHYS_SERN3D))
    # 重複ノード: 座標一致・ID 相異 (カウルスリット・側壁スリット)
    check("重複ノード数 = カウル + 側壁", info["nodes"] == info["ni"] * info["NJ"] * info["nz"] + info["n_dup_cowl"] + info["n_dup_side"])

# --- R4: ext_top (機体上面 + 自由流バンド) の押し出し ---------------------------------------------------
prm = SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, W=2.0, Z_ext=1.5, W_vehicle=3.0,
                       interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"], ext_top=True, nj_ext_top=9, top_depth=1.5,
                       vehicle_taper=0.3, vehicle_wedge_deg=3.0)
coords, hexes, B, info, y_mid = generate_sern_mesh3d(d, prm)
print(f"--- ext_top: cells {info['cells']} nodes {info['nodes']} top nodes {info['n_top_nodes']} vehicle_top faces {len(B['vehicle_top'])}")
nun, miss, extra = closure(hexes, B)
check("ext_top: 境界の閉性 (プルーム上線が内部面になり top_out は上面だけ)", miss == 0 and extra == 0, f"unshared {nun} missing {miss} extra {extra}")
ir = info["i_ramp_te"]
check("ext_top: vehicle_top 面数 = i_ramp_te × (nz−1)", len(B["vehicle_top"]) == ir * (info["nz"] - 1))
check("ext_top: top_out 面数 = (ni−1)(nz−1) (上面のみ)", len(B["top_out"]) == (info["ni"] - 1) * (info["nz"] - 1))
yt_faces = np.array([coords[list(q), 1].mean() for q in B["top_out"]]) / prm.scale
check("ext_top: top_out は y3 + top_depth より上", np.all(yt_faces > info["y_veh"] - 1e-9))
vt = np.array([coords[list(q), 1].mean() for q in B["vehicle_top"]]) / prm.scale
rp = np.array([coords[list(q), 1].mean() for q in B["ramp"]]) / prm.scale
check("ext_top: 機体上面はランプ (下面) より上", vt.min() >= rp.min() and vt.max() >= rp.max())
# 後縁でテーパが y_e に着地 (x = L_ramp の上面ノード = プルーム上線ノードと共有)
te_nodes = {n for q in B["vehicle_top"] for n in q if abs(coords[n, 0] / prm.scale - info["L_ramp"]) < 1e-9}
ramp_te = {n for g in ("ramp", "vehicle", "underside_far") for q in B[g] for n in q if abs(coords[n, 0] / prm.scale - info["L_ramp"]) < 1e-9}
check("ext_top: 後縁 station のノードは機体上面と下面 (ramp/vehicle/underside_far) で共有 (テーパが y_e に着地)", te_nodes and te_nodes == ramp_te, f"{len(te_nodes)} vs {len(ramp_te)}")
check("ext_top: hex は非退化 (体積 > 0)", np.all(np.abs(np.linalg.det(np.stack([coords[hexes[:, 1]] - coords[hexes[:, 0]], coords[hexes[:, 3]] - coords[hexes[:, 0]], coords[hexes[:, 4]] - coords[hexes[:, 0]]], axis=1))) > 1e-18))
try:
    generate_sern_mesh3d(d, SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, ext_top=True, vehicle_taper=0.0,
                                             interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"])); ok = False
except ValueError:
    ok = True
check("ext_top: vehicle_taper=0 (鉛直 base) は未実装として拒否", ok)

# --- ランプ丸め (ramp_fillet) の 3D 移植 ---
prm_f = SernMesh3DParams(ni_up=8, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, W=2.0, Z_ext=1.5, ramp_fillet=0.1,
                         interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"], ext_top=True, nj_ext_top=9, vehicle_taper=0.3)
coords, hexes, B, info, y_mid = generate_sern_mesh3d(d, prm_f)
nun, miss, extra = closure(hexes, B)
check("fillet: 境界の閉性", miss == 0 and extra == 0, f"missing {miss} extra {extra}")
xr = np.array([coords[list(q), 0].mean() for q in B["ramp"]]) / prm_f.scale; yr = np.array([coords[list(q), 1].mean() for q in B["ramp"]]) / prm_f.scale
near = (xr > -0.05) & (xr < 0.05)
check("fillet: 角部近傍のランプ面は y > 1 (流体と反対側へ膨らむ)", near.any() and np.all(yr[near] >= 1.0 - 1e-12), f"min y {yr[near].min():.4f}")
check("fillet: info に ramp_fillet", info["ramp_fillet"] == 0.1)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
