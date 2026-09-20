#!/usr/bin/env python3
"""⑤ SERN 3D hex メッシュ (mesh_sern3d) の整合テスト: 境界の閉性・タグ分割 (R2 の vehicle タグ)・重複ノード。"""
import sys
from collections import Counter
from pathlib import Path
import numpy as np
from dataclasses import replace
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
ir = info["i_ramp_te"]; k_sw = info["k_sw"]; nzf = info["nz"]
# R4e 案 (d): 機体上面は **幅内 (k < k_sw) かつ後縁より上流**だけ。幅外はバンド、後縁より下流は後流ブロックとの内部面
_exp = ir * k_sw
check("ext_top: vehicle_top 面数 = 後縁 station × 幅内 z セル数", len(B["vehicle_top"]) == _exp, f"{len(B['vehicle_top'])} vs {_exp}")
zc = np.array([coords[list(q), 2].mean() for q in B["vehicle_top"]]) / prm.scale
check("ext_top: 機体上面は全面が幅内 (z ≤ W/2)  ← 中央値でなく最大値で見る (codex plan-3 M5)",
      float(zc.max()) <= 0.5 * prm.W + 1e-9, f"max z {float(zc.max()):.3f}")
check("ext_top: 幅外の旧ランプ線は内部面 (vehicle / underside_far は 0)",
      len(B.get("vehicle", [])) + len(B.get("underside_far", [])) == 0,
      f"vehicle {len(B.get('vehicle',[]))} underside_far {len(B.get('underside_far',[]))}")
check("ext_top: 機体側面 (vehicle_side) が在る", len(B.get("vehicle_side", [])) > 0, f"{len(B.get('vehicle_side',[]))} faces")

# --- R4e の核心: 幅外を横断して塞ぐ壁が 1 面も無いこと ---
# 「幅外 (z > W/2) に在って法線が x 方向を向く壁面」= 旧実装が 224 面持っていた閉じ壁
def _face_normal_x(q):
    P = coords[list(q)]
    n = np.cross(P[1] - P[0], P[2] - P[0])
    nl = np.linalg.norm(n)
    return abs(n[0]) / nl if nl > 0 else 0.0
_wall_groups = ("vehicle_side", "vehicle_base", "vehicle_top", "ramp", "cowl_in", "cowl_out",
                "sidewall_in", "sidewall_out", "vehicle", "underside_far")
_blockers = [q for g in _wall_groups for q in B.get(g, [])
             if coords[list(q), 2].mean() / prm.scale > 0.5 * prm.W + 1e-9 and _face_normal_x(q) > 0.5]
check("ext_top: 幅外を横断する壁面が 0 (R4e: 旧実装は 224 面で 5.72 MPa を溜めた)",
      len(_blockers) == 0, f"{len(_blockers)} faces")

# --- 機体ベース: 幅内のみ・面積が t_base × W/2 と一致 ---
vb = B.get("vehicle_base", [])
check("ext_top: 機体ベース (vehicle_base) が在る", len(vb) > 0, f"{len(vb)} faces")
check("ext_top: ベースは幅内のみ (z ≤ W/2)",
      all(coords[list(q), 2].mean() / prm.scale <= 0.5 * prm.W + 1e-9 for q in vb))
def _quad_area(q):
    P = coords[list(q)]
    return 0.5 * (np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0])) + np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[0])))
A_base = sum(_quad_area(q) for q in vb) / prm.scale ** 2
A_exp = float(info["t_base"]) * 0.5 * prm.W
check("ext_top: ベース面積 = t_base × W/2 (設計値と照合)", abs(A_base - A_exp) < 1e-6 * max(A_exp, 1.0),
      f"{A_base:.6f} vs {A_exp:.6f}")
check("ext_top: top_out 面数 = (ni−1)(nz−1) (上面のみ)", len(B["top_out"]) == (info["ni"] - 1) * (info["nz"] - 1))
yt_faces = np.array([coords[list(q), 1].mean() for q in B["top_out"]]) / prm.scale
check("ext_top: top_out は y3 + top_depth より上", np.all(yt_faces > info["y_veh"] - 1e-9))
vt = np.array([coords[list(q), 1].mean() for q in B["vehicle_top"]]) / prm.scale
rp = np.array([coords[list(q), 1].mean() for q in B["ramp"]]) / prm.scale
check("ext_top: 機体上面はランプ (下面) より上", vt.min() >= rp.min() and vt.max() >= rp.max())

# --- 幾何: **符号付き** Jacobian (絶対値では負向き要素を見逃す, codex plan-3 M5) ---
_J = np.linalg.det(np.stack([coords[hexes[:, 1]] - coords[hexes[:, 0]],
                             coords[hexes[:, 3]] - coords[hexes[:, 0]],
                             coords[hexes[:, 4]] - coords[hexes[:, 0]]], axis=1))
check("ext_top: hex の符号付き Jacobian が全て同符号かつ非退化", np.all(_J > 1e-18) or np.all(_J < -1e-18),
      f"min {_J.min():.3e} max {_J.max():.3e} 負 {int((_J<0).sum())}/{len(_J)}")

# --- float32 変換後の節点衝突 (変換器は float32 で書く) ---
# スリット (カウル板厚 0 / 側壁) の重複ノードは**設計上わざと座標一致**しているので、
# 判定は「float32 にして **新たに** 衝突するノードが 0」かどうか
_u64 = np.unique(coords, axis=0).shape[0]
_u32 = np.unique(coords.astype(np.float32), axis=0).shape[0]
check("ext_top: float32 化で新たな節点衝突が起きない", _u32 == _u64, f"float32 {_u32} vs float64 {_u64}")
_nd = info["n_dup_cowl"] + info["n_dup_side"]
check("ext_top: 座標一致ノードはスリット由来だけ (カウル + 側壁)", coords.shape[0] - _u64 == _nd,
      f"{coords.shape[0] - _u64} vs {_nd}")

# --- 形状が格子に依存しないこと (codex plan-3 M3) ---
_p2 = replace(prm, nj_ext_top=13, first_top_frac=0.01)
_c2, _h2, _B2, _i2, _ = generate_sern_mesh3d(d, _p2)
check("ext_top: first_top_frac を変えても機体上面の高さが動かない",
      abs(_i2["y_veh"] - info["y_veh"]) < 1e-12, f"{_i2['y_veh']:.6f} vs {info['y_veh']:.6f}")
_n2, _m2, _e2 = closure(_h2, _B2)
check("ext_top: 細かい格子でも閉性", _m2 == 0 and _e2 == 0, f"missing {_m2} extra {_e2}")
try:
    generate_sern_mesh3d(d, replace(prm, first_top_frac=0.05)); ok = False
except ValueError:
    ok = True
check("ext_top: クリアランスに対して粗すぎる格子は生成を失敗させる (形状を動かさない)", ok)
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



# --- 壁第 1 層の x ブレンド (plan sern-3d §4.41) ---
_pb = SernMesh3DParams(ni_up=6, ni_noz=20, ni_plume=30, nj_top=15, nj_bot=11, nz_in=7, nz_out=6, W=2.0, Z_ext=1.5,
                       interface_angle=float(k.TH[-1, 0]), top_ext_angle=d.info["theta_e"], first_wall_frac=4.0e-4)
c0, h0, B0, i0, _ = generate_sern_mesh3d(d, _pb)
c1, h1, B1, i1, _ = generate_sern_mesh3d(d, replace(_pb, first_wall_frac_far=0.0))
check("xblend: 既定 (far=0) は座標がビット一致", np.array_equal(c0, c1))
c2, h2, B2, i2, _ = generate_sern_mesh3d(d, replace(_pb, first_wall_frac_far=4.0e-3, wall_frac_blend_len=0.5))
check("xblend: 節点数・要素数は変わらない", c2.shape == c0.shape and h2.shape == h0.shape)
check("xblend: 壁の形状は動かない (ランプ線 y の最大)",
      abs(float(np.max(c2[:, 1])) - float(np.max(c0[:, 1]))) < 1e-12,
      f"{np.max(c2[:,1]):.12f} vs {np.max(c0[:,1]):.12f}")
# 生成後の実座標で第 1 層厚を測る (入力値でなく結果を見る)
_Lc = float(d.cowl_xy[-1, 0]); _Lr = float(d.L_ramp)
def _first_layer(coords, xq):
    """x = xq の station で、上線 (ランプ) 直下の第一層厚を返す"""
    m = np.abs(coords[:, 0] - xq) < 1e-9
    ys = np.unique(np.round(coords[m, 1], 12))
    return float(ys[-1] - ys[-2])
_xs = np.unique(c0[:, 0])
_x_in = _xs[np.argmin(np.abs(_xs - 0.5 * _Lr))]              # 壁の内側
_x_far = _xs[np.argmin(np.abs(_xs - (_Lr + 1.5)))]           # ブレンド完了後
check("xblend: 壁の内側では第一層厚が first_wall_frac のまま",
      abs(_first_layer(c2, _x_in) - _first_layer(c0, _x_in)) < 1e-12,
      f"{_first_layer(c2,_x_in):.3e} vs {_first_layer(c0,_x_in):.3e}")
check("xblend: 壁の下流では第一層厚が粗くなる (>= 5 倍)",
      _first_layer(c2, _x_far) >= 5.0 * _first_layer(c0, _x_far),
      f"{_first_layer(c2,_x_far):.3e} vs {_first_layer(c0,_x_far):.3e}")
_n2, _m2, _e2 = closure(h2, B2)
check("xblend: 境界の閉性", _m2 == 0 and _e2 == 0, f"missing {_m2} extra {_e2}")
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
