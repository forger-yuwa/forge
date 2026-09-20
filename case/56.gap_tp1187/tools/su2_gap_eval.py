#!/usr/bin/env python3
"""SU2 の VTU からすきま壁の熱流束を深さ別に取り出す (forge との重ね合わせ用)。

`procedures/su2-cross-check.md` の「ライン比較プロトコル」に従い、**同一メッシュ・同一条件**の
forge 結果と同じ座標で比べる。SU2 の `Heat_Flux` は壁面のみ有限で、内部節点では 0。
"""
import argparse, sys
from pathlib import Path
import re
import numpy as np

CASE = Path(__file__).resolve().parents[1]

_DT = {"Float32": np.float32, "Float64": np.float64,
       "Int32": np.int32, "Int64": np.int64, "UInt8": np.uint8,
       "UInt32": np.uint32, "UInt64": np.uint64}


def read_vtu_appended(path):
    """SU2 の appended-raw VTU を自前で読む。

    meshio は `buffer size must be a multiple of element size` で落ち、vtk の XML パーサも
    生バイナリ部で `not well-formed` になる (2026-09-20 実測)。形式は単純で、
    `<AppendedData encoding="raw">_` の後に [header_type の byte 数][データ] が並ぶだけ。
    """
    raw = Path(path).read_bytes()
    i = raw.index(b"<AppendedData")
    head = raw[:i].decode("latin1")
    j = raw.index(b"_", i) + 1
    blob = raw[j:raw.index(b"</AppendedData>", j)]
    ht = _DT[re.search(r'header_type="(\w+)"', head).group(1)]
    hsz = np.dtype(ht).itemsize

    def get(off, dtype, ncomp):
        n = int(np.frombuffer(blob[off:off + hsz], dtype=ht, count=1)[0])
        v = np.frombuffer(blob[off + hsz:off + hsz + n], dtype=dtype)
        return v.reshape(-1, ncomp) if ncomp > 1 else v

    arrays, points = {}, None
    for m in re.finditer(r"<DataArray([^>]*)/>", head):
        at = m.group(1)
        nm = (re.search(r'Name="([^"]*)"', at) or [None, ""])[1]
        dt = _DT[re.search(r'type="(\w+)"', at).group(1)]
        nc = int((re.search(r'NumberOfComponents=\s*"(\d+)"', at) or [None, "1"])[1])
        off = int(re.search(r'offset="(\d+)"', at).group(1))
        v = get(off, dt, nc)
        if nm == "" and points is None:
            points = v
        elif nm:
            arrays[nm] = v
    return points, arrays



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtu")
    ap.add_argument("--w", type=float, default=1.8e-3)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    pts, pd = read_vtu_appended(a.vtu)
    c = pts[:, :2]
    q = np.asarray(pd["Heat_Flux"]).ravel()
    T = np.asarray(pd["Temperature"]).ravel()
    xd = 0.9e-3
    # 下流側 (前向き壁) の壁節点: x = +W/2 (鉛直部) と下流円弧
    wall = np.abs(q) > 0
    sel = wall & (c[:, 0] > 1e-9) & (c[:, 1] < 1e-9)
    d = -c[sel, 1]; qq = q[sel]; o = np.argsort(d)
    d, qq = d[o], qq[o]
    print(f"{Path(a.vtu).name}: 壁節点 {int(wall.sum())}, 下流側すきま壁 {len(d)}")
    print(f"{'z/W':>6} {'深さ[mm]':>9} {'q_w [W/m2]':>13} {'T_w [K]':>9}")
    for z in (1.39, 2.83, 4.22, 8.44, 14.11, 21.17):
        if z * a.w > d.max():
            continue
        print(f"{z:6.2f} {z*a.w*1e3:9.2f} {np.interp(z*a.w, d, qq):13.4e} "
              f"{np.interp(z*a.w, d, T[sel][o]):9.3f}")
    if a.out:
        np.savetxt(a.out, np.column_stack([d, qq]), delimiter=",",
                   header="depth_m,q_w_W_m2", comments="")
        print(f"  → {a.out}")
    # 深部の温度分解能 (SU2 は倍精度)
    dep = -c[:, 1]
    print(f"\n{'z/W 帯':>10} {'T 幅 [K]':>12} {'|U| 最大 [m/s]':>15}")
    V = np.asarray(pd["Velocity"])[:, :2]
    for lo, hi in ((1, 3), (3, 5), (5, 10), (10, 35)):
        mm = (dep > lo * a.w) & (dep <= hi * a.w) & (np.abs(c[:, 0]) <= 0.9e-3 + 1e-9)
        if mm.sum():
            print(f"{lo:4d}-{hi:<5d} {float(T[mm].max()-T[mm].min()):12.3e} "
                  f"{np.hypot(V[mm,0],V[mm,1]).max():15.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
