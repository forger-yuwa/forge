#!/usr/bin/env python3
"""対称面 (既定 z = 0) で形状と値を鏡像反転し、全スパンの場を復元する。

半スパンで解いた 3D の結果を可視化するときに使う。ParaView の Reflect フィルタと違い、
**ベクトル成分の符号を正しく反転する**のが要点。

反転規則 (z → −z):
  スカラー (ro, P, T, k, omega, wall_dist ...)          そのまま
  ベクトルの z 成分 (Uz, roUz ...)                       符号反転
  スカラーの z 微分 (dPdz, drodz, dKdz, dOmegadz ...)    符号反転
  ベクトルの微分 (dUidxj)  z が**奇数回**現れる成分だけ反転
      dUxdz, dUydz, dUzdx, dUzdy → 反転    dUzdz → そのまま

対称面上のノードは重複させない (座標一致ノードは最近傍補間や node 方式で事故の元になる)。

使い方:
  python3 mirror_symmetry.py res_12000.h5 -o res_12000_full.h5      # 体積・面のどちらでも可
  python3 mirror_symmetry.py res_ramp_4_500.h5 --axis y --plane 0.0
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import numpy as np

AXIS = {"x": 0, "y": 1, "z": 2}
_CONNE_LEN = 0


# --- 変数の型 (codex plan レビュー 2 の M4/M5 採用) --------------------------------------------
# 名前の末尾で判定すると `limiter_Uz` (制限係数=スカラー) を速度成分と誤判定する。
# 反射 R = diag(1,1,-1) に対し、極性ベクトルは Rv、速度勾配は R G Rᵀ、軸性ベクトル (渦度) は det(R)·Rv。
POLAR = {"Ux", "Uy", "Uz", "roUx", "roUy", "roUz", "twall_x", "twall_y", "twall_z"}   # 極性ベクトル成分
AXIAL_PREFIX = ("vort", "omega_vec")        # 軸性ベクトル: 反射で x,y が反転し z は維持
SCALAR_EXACT = {"ro", "P", "T", "e", "h0", "k", "omega", "roK", "roOmega", "roe", "sonic",
                "vis_lam", "vis_turb", "wall_dist", "volume", "dt_local", "cfl", "patch",
                "delta_les", "l_des", "rd_des", "fd_shield", "fe_iddes", "ducros", "Taw_diag",
                "thermCond", "Pk_diag", "wf_pk", "roK_wf", "wf_irep_flag", "axisym_divU", "condTinvFail",
                # 壁面・境界出力のスカラー (面出力 res_<tag>_<pid>_<step>.h5)
                "Ps", "Pt", "Ts", "Tt", "kb", "omegab", "qwall", "utau", "ypls"}
SCALAR_PREFIX = ("limiter_", "src_jac_", "transport_diag_", "res_", "Y", "roY", "X", "roX")


def _axis_count(comp: str, ax: str) -> int:
    return 1 if comp == ax else 0


def classify(name: str, ax: str):
    """(反転するか, 分類) を返す。分類が 'unknown' なら呼び出し側で警告する。"""
    # 勾配 d<var>d<xyz>
    m = re.fullmatch(r"d(.+)d([xyz])", name)
    if m:
        var, wrt = m.group(1), m.group(2)
        n = _axis_count(wrt, ax)
        if var in POLAR or (len(var) >= 2 and var[-1] in AXIS and var[:-1] + "x" in POLAR):
            n += _axis_count(var[-1], ax)
            return n % 2 == 1, "grad_vector"
        return n % 2 == 1, "grad_scalar"     # スカラーの勾配
    if name in POLAR:
        return name[-1] == ax, "polar"
    if name.startswith(AXIAL_PREFIX):
        return name[-1] != ax, "axial"       # 軸性: 法線成分は維持、面内成分が反転
    if name in SCALAR_EXACT or name.startswith(SCALAR_PREFIX):
        return False, "scalar"
    return False, "unknown"


def flip_sign(name: str, ax: str) -> bool:
    return classify(name, ax)[0]


# CONNE の先頭は**節点数ではなく XDMF の要素コード** (res_h5_to_vtu.py と同じ表)
XDMF_NNODE = {4: 3, 5: 4, 6: 4, 7: 5, 8: 6, 9: 8}   # tri / quad / tetra / pyramid / wedge / hex
# 鏡像は向きを反転するので、要素種別ごとに**向きを戻す頂点置換**を当てる。
# `ids[::-1]` は quad では正しいが hex では Jacobian が負のままになる (codex M5)。
MIRROR_PERM = {
    4: (0, 2, 1),                               # tri: 2 点入れ替え
    5: (0, 3, 2, 1),                            # quad: 逆回り
    6: (0, 1, 3, 2),                            # tetra: 2 点入れ替え
    7: (0, 3, 2, 1, 4),                         # pyramid: 底面を逆回り
    8: (3, 4, 5, 0, 1, 2),                      # wedge: 2 つの三角形を入れ替え
    9: (4, 5, 6, 7, 0, 1, 2, 3),                # hex: 底面と上面を入れ替え
}


def _read(src: Path):
    with h5py.File(src, "r") as f:
        return (np.asarray(f["MESH"]["COORD"][:]).reshape(-1, 3),
                np.asarray(f["MESH"]["CONNE"][:]),
                {k: np.asarray(f["VALUE"][k][:]) for k in f["VALUE"].keys()})


def merge(srcs) -> tuple:
    """複数の面出力を 1 つのメッシュに束ねる。`patch` に面ごとの整数 id を入れる
    (ParaView で色分けすれば、どの壁がどこかが一目で分かる)。共通しない変数は NaN 埋め。"""
    coords, connes, vals, patch, off = [], [], [], [], 0
    keys = set()
    for i, sp in enumerate(srcs):
        c, cn, v = _read(Path(sp))
        # CONNE のノード id をオフセット
        out, j = [], 0
        while j < len(cn):
            code = int(cn[j]); n = XDMF_NNODE[code]
            out.append(code); out.extend(int(x) + off for x in cn[j + 1:j + 1 + n]); j += n + 1
        coords.append(c); connes.append(np.asarray(out, dtype=cn.dtype)); vals.append(v)
        patch.append(np.full(len(c), i, dtype=np.float32)); off += len(c); keys |= set(v)
    val = {}
    for k in sorted(keys):
        val[k] = np.concatenate([v[k] if k in v and len(v[k]) == len(c) else np.full(len(c), np.nan, np.float32)
                                 for v, c in zip(vals, coords)])
    val["patch"] = np.concatenate(patch)
    return np.vstack(coords), np.concatenate(connes), val


def mirror(src, dst: Path, ax: str = "z", plane: float = 0.0, tol: float = 1e-9) -> dict:
    a = AXIS[ax]
    if isinstance(src, (list, tuple)):
        coord, conne, val = merge(src)
    else:
        coord, conne, val = _read(Path(src))

    nN = len(coord)
    on = np.abs(coord[:, a] - plane) <= tol          # 対称面上 = 重複させない
    keep = ~on
    nNew = int(keep.sum())
    # 旧 id -> 鏡像側 id。対称面上は自分自身を指す
    newid = np.full(nN, -1, dtype=np.int64)
    newid[keep] = nN + np.arange(nNew)
    newid[on] = np.where(on)[0]

    mcoord = coord[keep].copy()
    mcoord[:, a] = 2.0 * plane - mcoord[:, a]
    coord_out = np.vstack([coord, mcoord])

    # CONNE: 要素ごとにノード id を写し、**向きを保つため頂点順を反転**する
    out, i, ncell = [], 0, 0
    while i < len(conne):
        code = int(conne[i])
        if code not in XDMF_NNODE:
            raise ValueError(f"未対応の XDMF 要素コード {code} (CONNE の先頭は節点数ではなくコード)")
        n = XDMF_NNODE[code]; ids = conne[i + 1:i + 1 + n]
        out.append(code); out.extend(int(newid[ids[q]]) for q in MIRROR_PERM[code])
        i += n + 1; ncell += 1
    conne_out = np.concatenate([conne, np.asarray(out, dtype=conne.dtype)])

    flipped, unknown, skipped = [], [], []
    val_out = {}
    for k, v in val.items():
        if len(v) != nN:                              # ノード数と合わない配列 (セル中心など) は型が判らないので複製しない
            skipped.append(k); continue
        fl, kind = classify(k, ax)
        if kind == "unknown":
            unknown.append(k)
        mv = v[keep].copy()
        if fl:
            mv = -mv; flipped.append(k)
        val_out[k] = np.concatenate([v, mv])

    with h5py.File(src if not isinstance(src, (list, tuple)) else src[0], "r") as f0:
        root_attrs = dict(f0.attrs)
        val_attrs = {k: dict(f0["VALUE"][k].attrs) for k in f0["VALUE"].keys() if f0["VALUE"][k].attrs}
    with h5py.File(dst, "w") as g:
        for ak, av in root_attrs.items():
            g.attrs[ak] = av                          # h0_includes_k 等を落とさない
        mg = g.create_group("MESH")
        mg.create_dataset("COORD", data=coord_out.reshape(-1).astype(coord.dtype))
        mg.create_dataset("CONNE", data=conne_out)
        vg = g.create_group("VALUE")
        for k, v in val_out.items():
            d = vg.create_dataset(k, data=v)
            for ak, av in val_attrs.get(k, {}).items():
                d.attrs[ak] = av
    global _CONNE_LEN
    _CONNE_LEN = len(conne_out)
    return {"nodes": nN, "nodes_out": len(coord_out), "on_plane": int(on.sum()),
            "cells": ncell, "cells_out": ncell * 2, "flipped": sorted(flipped),
            "unknown": sorted(unknown), "skipped": sorted(skipped), "keys": sorted(val_out)}


def write_xmf(h5: Path, nnode: int, ncell: int, keys, time: float = 0.0) -> Path:
    """対の .xmf を書く (ParaView は .h5 単体では開けない)。元の res_*.xmf と同じ書式。"""
    xmf = h5.with_suffix(".xmf")
    a = "".join(f"        <Attribute Name='{k}' Center='Node' >\n"
                f"          <DataItem Format='HDF' DataType='Float' Dimensions='{nnode}'>\n"
                f"            {h5.name}:VALUE/{k}\n          </DataItem>\n        </Attribute>\n" for k in keys)
    xmf.write_text(
        "<?xml version='1.0' ?>\n<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n<Xdmf>\n  <Domain>\n"
        "    <Grid GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n"
        f"    <Time TimeType='Single' Value='{time:g}'/>\n      <Grid Name='mirrored'>\n"
        f"        <Topology Type='Mixed' NumberOfElements='{ncell}'>\n"
        f"          <DataItem Format='HDF' DataType='Int' Dimensions='{_CONNE_LEN}'>\n"
        f"            {h5.name}:MESH/CONNE\n          </DataItem>\n        </Topology>\n"
        "        <Geometry Type='XYZ'>\n"
        f"          <DataItem Format='HDF' DataType='Float' Dimensions='{nnode * 3}'>\n"
        f"            {h5.name}:MESH/COORD\n          </DataItem>\n        </Geometry>\n"
        + a + "      </Grid>\n    </Grid>\n  </Domain>\n</Xdmf>\n")
    return xmf


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="対称面で形状と値を鏡像反転して全スパンを復元する")
    ap.add_argument("src", nargs="+", help="res_*.h5 (複数指定するとまとめて 1 つにする)")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--axis", default="z", choices=list(AXIS))
    ap.add_argument("--plane", type=float, default=0.0)
    ap.add_argument("--tol", type=float, default=1e-9)
    a = ap.parse_args(argv)
    if len(a.src) == 1:
        src = Path(a.src[0]); dst = Path(a.out) if a.out else src.with_name(src.stem + "_full" + src.suffix)
        label = src.name
    else:
        src = [Path(x) for x in a.src]; dst = Path(a.out) if a.out else src[0].with_name("walls_full.h5")
        label = f"{len(src)} 面"
        print("  まとめる面 (patch id):")
        for i, sp in enumerate(src):
            print(f"    {i}  {Path(sp).name}")
    info = mirror(src, dst, a.axis, a.plane, a.tol)
    xmf = write_xmf(dst, info["nodes_out"], info["cells_out"], info["keys"])
    print(f"{label} -> {dst.name} (+ {xmf.name})")
    print(f"  ノード {info['nodes']} -> {info['nodes_out']}  (対称面上 {info['on_plane']} は複製せず)")
    print(f"  要素   {info['cells']} -> {info['cells_out']}")
    print(f"  符号反転した変数 ({len(info['flipped'])}): {', '.join(info['flipped']) or 'なし'}")
    if info["unknown"]:
        print(f"  ★ 型が未登録 (そのまま複製した。符号規則の確認が要る): {', '.join(info['unknown'])}")
    if info["skipped"]:
        print(f"  ノード数と合わず除外: {', '.join(info['skipped'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
