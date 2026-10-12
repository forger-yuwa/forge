r"""解析メッシュの構造ファイル `mesh_structure.json` (plans/active/tooling-nozzle-wt-deliverables.md §4.3・§6 V3)。

run の `nozzle.msh` (`meshing/mesh2d.py` の `write_msh41_2d` が書いた msh4.1) を読み、STEP の壁曲線から gmsh で
**同じ構造**のメッシュを作り直すための情報 (station $x_i$・半径方向の比 $s_{ij}$・接続の規則・物理名とタグ) を取り出す。
メッシャは再実装しない (メッシャはブランチごとに違い、今後も変わるため)。作り直しは単体スクリプト
`rebuild_from_step.py` (gmsh、forge_design に依存しない) が行う。

- 接続の規則 `mesh2d_ij_v1`: 節点番号 $(i,j)\to i\,n_j+j$ ($i$ = station、軸方向 / $j$ = 軸 0 → 壁 $n_j-1$)、
  四角形は $i$ が外側・$j$ が内側のループで $(i,j),(i{+}1,j),(i{+}1,j{+}1),(i,j{+}1)$ (反時計回り)、
  境界辺は軸 $j=0$・壁 $j=n_j-1$ ($i$ 増加の向き)、入口 $i=0$・出口 $i=n_i-1$ ($j$ 増加の向き)。
  要素番号は四角形が先で、続いて境界辺のブロックを msh に現れた順 (`write_msh41_2d` は inlet・outlet・wall・axis)。
- 縦線の検査: 各 station で $\max_j |x_{ij}-x_{i0}| \le 10^{-9} L$ ($L$ = 全長)。外れたら作らない (壁法線の近壁層など)。
- 比 $s_{ij} = r_{ij}/r_{i,n_j-1}$。全 station の行が 1e-15 以内で同じなら `s_common` (1 行)、違えば `s` ($n_i\times n_j$)。
  msh は座標を 10 有効桁で書くので、実際の msh では行が丸めの分だけ違い、ふつうは `s` になる。軸の比 $s_{i0}$ は 0 に置く
  (軸の節点が $r=0$ であることを検査したうえで)。
- 単位: msh の座標は m、構造ファイルの長さは mm (STEP と同じ)。

usage:
  design/.venv-opt/bin/python -m forge_design.export.mesh_structure NOZZLE_MSH [--out mesh_structure.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

FORMAT = "forge.mesh_structure/v1"
CONNECTIVITY = "mesh2d_ij_v1"
VERTICAL_TOL_REL = 1e-9      # 縦線の許容 (全長に対する比; plan §4.3)
S_COMMON_TOL = 1e-15         # 比の行が全 station で同じとみなす許容
AXIS_TOL_REL = 1e-12         # 軸の節点の |r| と全節点の |z| の許容 (全長に対する比)
M_TO_MM = 1000.0
SIDES = ("inlet", "outlet", "wall", "axis")   # 境界辺の名前 (msh の物理名) = 辺の位置


class StructureError(ValueError):
    """msh が `mesh2d_ij_v1` の構造 (縦線・節点番号の規則・軸) に従わない。メッセージに理由を書く。"""


# --------------------------------------------------------------------------------------------
# msh4.1 の読み込み (write_msh41_2d の書式だけ)
# --------------------------------------------------------------------------------------------

def _section(lines: list[str], name: str) -> list[str]:
    """`$name` 〜 `$Endname` の中身の行。無い・重複は ValueError。"""
    beg = [k for k, s in enumerate(lines) if s.strip() == f"${name}"]
    end = [k for k, s in enumerate(lines) if s.strip() == f"$End{name}"]
    if len(beg) != 1 or len(end) != 1 or end[0] <= beg[0]:
        raise ValueError(f"msh: 節 ${name} が 1 つでない (開始 {len(beg)}・終了 {len(end)})")
    return lines[beg[0] + 1:end[0]]


def read_msh41_mesh2d(path):
    """`write_msh41_2d` の書式の msh4.1 (ASCII) を読む純 Python の読み込み器。

    書式の条件 (違えば ValueError): `$MeshFormat 4.1 0 8`・物理名 (曲線 4 + 面 1)・実体は点 0・曲線 4・面 1・体 0 で、
    それぞれ物理タグ 1 つ・節点ブロックが 1 つ (面に載る、助変数なし、番号は 1..N の並べ替え)・要素ブロックは
    最初が面の四角形 (型 3)、続いて曲線ごとの線分 (型 1) が 4 つ、要素番号はファイルの順に 1..M。

    戻り値: (coords (N,3) float64 [m], quads (M,4) int64 0 始まり, 境界辺 dict 名前 → (K,2) int64 0 始まり
    (msh に現れた順), 物理名 → タグ dict)。境界辺の名前は曲線の物理名。
    """
    text = Path(path).read_text(encoding="utf-8")
    lines = text.splitlines()
    fmt = _section(lines, "MeshFormat")
    if len(fmt) < 1 or fmt[0].split() != ["4.1", "0", "8"]:
        raise ValueError(f"msh: 書式が 4.1 ASCII でない ({fmt[:1]})")
    # 物理名
    pn = _section(lines, "PhysicalNames")
    npn = int(pn[0])
    if npn != len(pn) - 1:
        raise ValueError("msh: $PhysicalNames の件数が合わない")
    phys: dict[str, int] = {}
    phys_dim: dict[int, tuple[int, str]] = {}
    for s in pn[1:]:
        a = s.split(maxsplit=2)
        dim, tag, name = int(a[0]), int(a[1]), a[2].strip().strip('"')
        if name in phys or tag in phys_dim:
            raise ValueError(f"msh: 物理名・タグが重複 ({name!r}, {tag})")
        phys[name] = tag
        phys_dim[tag] = (dim, name)
    # 実体
    ent = _section(lines, "Entities")
    cnt = [int(v) for v in ent[0].split()]
    if cnt != [0, 4, 1, 0]:
        raise ValueError(f"msh: 実体の数が 点 0・曲線 4・面 1・体 0 でない ({cnt})")
    curve_phys: dict[int, int] = {}
    for s in ent[1:5]:
        a = s.split()
        if int(a[7]) != 1:
            raise ValueError(f"msh: 曲線 {a[0]} の物理タグが 1 つでない")
        curve_phys[int(a[0])] = int(a[8])
    a = ent[5].split()
    if int(a[7]) != 1:
        raise ValueError("msh: 面の物理タグが 1 つでない")
    surf_tag, surf_phys = int(a[0]), int(a[8])
    for ptag, dim in [(t, 1) for t in curve_phys.values()] + [(surf_phys, 2)]:
        if ptag not in phys_dim or phys_dim[ptag][0] != dim:
            raise ValueError(f"msh: 物理タグ {ptag} が $PhysicalNames に無いか次元 ({dim}) が違う")
    if len(set(curve_phys.values())) != 4:
        raise ValueError("msh: 4 本の曲線の物理タグが重複")
    # 節点
    nd = _section(lines, "Nodes")
    nblk, nnode, tmin, tmax = (int(v) for v in nd[0].split())
    if nblk != 1:
        raise ValueError(f"msh: 節点ブロックが 1 つでない ({nblk})")
    bd, bt, bpar, bn = (int(v) for v in nd[1].split())
    if (bd, bt, bpar, bn) != (2, surf_tag, 0, nnode):
        raise ValueError(f"msh: 節点ブロックが面 {surf_tag} の助変数なしの全節点でない ({nd[1]!r})")
    if len(nd) != 2 + 2 * nnode:
        raise ValueError("msh: 節点の行数が合わない")
    tags = np.array(nd[2:2 + nnode], dtype=np.int64)
    xyz = np.array(" ".join(nd[2 + nnode:]).split(), dtype=np.float64)
    if xyz.size != 3 * nnode:
        raise ValueError("msh: 節点座標が 3 成分でない")
    if (tmin, tmax) != (1, nnode) or not np.array_equal(np.sort(tags), np.arange(1, nnode + 1)):
        raise ValueError("msh: 節点番号が 1..N の並べ替えでない")
    coords = np.empty((nnode, 3))
    coords[tags - 1] = xyz.reshape(nnode, 3)
    # 要素
    el = _section(lines, "Elements")
    nblk, nel, emin, emax = (int(v) for v in el[0].split())
    if nblk != 5:
        raise ValueError(f"msh: 要素ブロックが 5 つ (四角形 1 + 線分 4) でない ({nblk})")
    k = 1
    quads = None
    bedges: dict[str, np.ndarray] = {}
    next_tag = 1
    for b in range(nblk):
        ed, et, etype, en = (int(v) for v in el[k].split())
        rows = el[k + 1:k + 1 + en]
        k += 1 + en
        npe = {3: 4, 1: 2}.get(etype)
        if b == 0 and (ed, et, etype) != (2, surf_tag, 3):
            raise ValueError("msh: 最初の要素ブロックが面の四角形 (型 3) でない")
        if b > 0 and (ed != 1 or etype != 1 or et not in curve_phys):
            raise ValueError(f"msh: 要素ブロック {b} が曲線の線分 (型 1) でない")
        arr = np.array(" ".join(rows).split(), dtype=np.int64)
        if len(rows) != en or arr.size != en * (1 + npe):
            raise ValueError(f"msh: 要素ブロック {b} の行数・節点数が合わない")
        arr = arr.reshape(en, 1 + npe)
        if not np.array_equal(arr[:, 0], np.arange(next_tag, next_tag + en)):
            raise ValueError("msh: 要素番号がファイルの順に 1..M でない")
        next_tag += en
        conn = arr[:, 1:] - 1
        if conn.min(initial=0) < 0 or conn.max(initial=0) >= nnode:
            raise ValueError(f"msh: 要素ブロック {b} が範囲外の節点を参照")
        if b == 0:
            quads = conn
        else:
            name = phys_dim[curve_phys[et]][1]
            if name in bedges:
                raise ValueError(f"msh: 曲線の要素ブロックが重複 ({name})")
            bedges[name] = conn
    if k != len(el) or (emin, emax, nel) != (1, next_tag - 1, next_tag - 1):
        raise ValueError("msh: $Elements の件数が合わない")
    return coords, quads, bedges, phys


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------------
# 構造の取り出し
# --------------------------------------------------------------------------------------------

def expected_connectivity(ni: int, nj: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """規則 `mesh2d_ij_v1` の四角形と境界辺 (`generate_axisym_mesh` と同じ並び)。"""
    nid = np.arange(ni * nj, dtype=np.int64).reshape(ni, nj)
    quads = np.stack([nid[:-1, :-1], nid[1:, :-1], nid[1:, 1:], nid[:-1, 1:]], axis=-1).reshape(-1, 4)
    bedges = {
        "inlet": np.stack([nid[0, :-1], nid[0, 1:]], axis=-1),
        "outlet": np.stack([nid[-1, :-1], nid[-1, 1:]], axis=-1),
        "wall": np.stack([nid[:-1, -1], nid[1:, -1]], axis=-1),
        "axis": np.stack([nid[:-1, 0], nid[1:, 0]], axis=-1),
    }
    return quads, bedges


def _first_mismatch(a: np.ndarray, b: np.ndarray) -> str:
    if a.shape != b.shape:
        return f"形 {a.shape} ≠ 規則 {b.shape}"
    k = int(np.argmax(np.any(a != b, axis=1)))
    return f"{k} 番目 (0 始まり) が {a[k].tolist()}、規則では {b[k].tolist()}"


def extract_structure(coords, quads, bedges, ni=None, nj=None, *, phys=None, msh_sha256=None,
                      source=None) -> dict:
    """msh の配列 → 構造 dict (`forge.mesh_structure/v1`)。規則に従わなければ StructureError (理由つき)。

    $n_i$ = 壁辺数 + 1、$n_j$ = 入口辺数 + 1 で推定する (引数で与えたときは一致を確かめる)。
    検査: (a) 節点番号・四角形・境界辺が規則 `mesh2d_ij_v1` と完全一致、(b) 縦線、(c) 軸の節点が $r=0$
    (と全節点が $z=0$)、ほかに station が $x$ 方向に狭義増加・各 station で $r$ が $j$ について狭義増加。
    phys (物理名 → タグ) を省くと `mesh2d.PHYS` の規約を使う。
    """
    coords = np.asarray(coords, dtype=np.float64)
    quads = np.asarray(quads, dtype=np.int64)
    bed = {k: np.asarray(v, dtype=np.int64).reshape(-1, 2) for k, v in bedges.items()}
    if set(bed) != set(SIDES):
        raise StructureError(f"境界辺の名前が {SIDES} でない ({sorted(bed)})")
    ni_e, nj_e = len(bed["wall"]) + 1, len(bed["inlet"]) + 1
    if ni is not None and int(ni) != ni_e:
        raise StructureError(f"ni = {ni} が壁辺数 + 1 = {ni_e} と違う")
    if nj is not None and int(nj) != nj_e:
        raise StructureError(f"nj = {nj} が入口辺数 + 1 = {nj_e} と違う")
    ni, nj = ni_e, nj_e
    if ni < 2 or nj < 3:
        raise StructureError(f"ni = {ni}・nj = {nj} が小さすぎる (ni ≥ 2・nj ≥ 3)")
    if coords.shape != (ni * nj, 3):
        raise StructureError(f"節点数 {coords.shape[0]} が ni·nj = {ni * nj} でない")
    # (a) 接続
    q_exp, b_exp = expected_connectivity(ni, nj)
    if not np.array_equal(quads, q_exp):
        raise StructureError("接続の検査: 四角形の節点番号が規則 (i,j)→i·nj+j と違う: " + _first_mismatch(quads, q_exp))
    for nm in SIDES:
        if not np.array_equal(bed[nm], b_exp[nm]):
            raise StructureError(f"接続の検査: 境界辺 {nm} が規則と違う: " + _first_mismatch(bed[nm], b_exp[nm]))
    X = coords[:, 0].reshape(ni, nj)
    R = coords[:, 1].reshape(ni, nj)
    Z = coords[:, 2]
    xi = X[:, 0].copy()
    L = float(xi[-1] - xi[0])
    if not (np.isfinite(coords).all() and L > 0):
        raise StructureError(f"座標に非有限値があるか全長 L = {L} ≤ 0")
    # (b) 縦線
    dev = np.abs(X - xi[:, None]).max(axis=1)
    vmax = float(dev.max()) / L
    if vmax > VERTICAL_TOL_REL:
        i = int(np.argmax(dev))
        j = int(np.argmax(np.abs(X[i] - xi[i])))
        raise StructureError(
            f"縦線の検査: station {i} (x = {xi[i] * M_TO_MM:.9g} mm) で節点 j = {j} の x が {dev[i] * M_TO_MM:.3e} mm "
            f"(= {vmax:.3e} L) ずれる (許容 {VERTICAL_TOL_REL:g} L)。格子線が縦でない (壁法線の近壁層など) 格子は v1 の構造ファイルで表せない")
    if not np.all(np.diff(xi) > 0):
        i = int(np.argmax(~(np.diff(xi) > 0)))
        raise StructureError(f"station の x が狭義増加でない (i = {i} → {i + 1})")
    # (c) 軸と平面
    amax = float(np.abs(R[:, 0]).max())
    if amax > AXIS_TOL_REL * L:
        i = int(np.argmax(np.abs(R[:, 0])))
        raise StructureError(f"軸の検査: 軸の節点 (j = 0) が r = 0 でない (station {i} で |r| = {amax:.3e} m、許容 {AXIS_TOL_REL:g} L)")
    zmax = float(np.abs(Z).max())
    if zmax > AXIS_TOL_REL * L:
        raise StructureError(f"平面の検査: 節点の z が 0 でない (|z| 最大 {zmax:.3e} m)")
    rw = R[:, -1].copy()
    if not np.all(rw > 0):
        raise StructureError("壁の半径が 0 以下の station がある")
    if not np.all(np.diff(R, axis=1) > 0):
        i = int(np.argmax(~np.all(np.diff(R, axis=1) > 0, axis=1)))
        raise StructureError(f"station {i} で r が j について狭義増加でない (格子が折れている)")
    S = R / rw[:, None]
    S[:, 0] = 0.0
    S[:, -1] = 1.0
    if phys is None:
        from forge_design.meshing.mesh2d import PHYS
        phys = dict(PHYS)
    surf = [k for k in phys if k not in SIDES]
    if not set(SIDES) <= set(phys) or len(surf) != 1:
        raise StructureError(f"物理名が {SIDES} + 面 1 つでない ({sorted(phys)})")
    physical = {nm: {"dim": 1, "tag": int(phys[nm])} for nm in SIDES}
    physical[surf[0]] = {"dim": 2, "tag": int(phys[surf[0]])}
    st: dict = {
        "format": FORMAT,
        "connectivity": CONNECTIVITY,
        "connectivity_note": "節点 (i,j) → i·nj + j (0 始まり; msh の番号は +1)。四角形 (i,j),(i+1,j),(i+1,j+1),(i,j+1) を i 外・j 内の順。"
                             "境界辺 axis j=0・wall j=nj-1 (i 増加)・inlet i=0・outlet i=ni-1 (j 増加)",
        "ni": int(ni),
        "nj": int(nj),
        "units": {"length": "mm", "msh_coords": "m", "mm_per_msh_unit": M_TO_MM, "s": "1 (r / r_wall)"},
        "physical": physical,
        "element_tag_order": [surf[0]] + [nm for nm in bedges],
        "length_mm": L * M_TO_MM,
        "x_mm": (xi * M_TO_MM).tolist(),
    }
    if float(np.abs(S - S[0][None, :]).max()) <= S_COMMON_TOL:
        st["s_common"] = S[0].tolist()
    else:
        st["s"] = S.tolist()
    # 作り直しの基準値 (元の msh の値)。比較の許容は第 1 層厚の 1 % (plan §6 V3)
    st["wall_r_mm"] = (rw * M_TO_MM).tolist()
    st["first_layer_mm"] = ((R[:, -1] - R[:, -2]) * M_TO_MM).tolist()
    st["reference_note"] = "wall_r_mm・first_layer_mm は元の msh の各 station の壁半径と第 1 層厚 (作り直しの基準値)"
    st["checks"] = {"vertical_max_dev_over_L": vmax, "vertical_tol_over_L": VERTICAL_TOL_REL,
                    "axis_max_abs_r_m": amax, "max_abs_z_m": zmax}
    st["source"] = {"msh": None if source is None else str(source), "sha256": msh_sha256}
    return st


def structure_from_msh(path) -> dict:
    """msh を読んで構造 dict を作る (sha256・物理名つき)。"""
    coords, quads, bedges, phys = read_msh41_mesh2d(path)
    return extract_structure(coords, quads, bedges, phys=phys, msh_sha256=sha256_file(path), source=Path(path).name)


def structure_s(struct: dict) -> np.ndarray:
    """構造 dict の比 (n_i, n_j)。"""
    ni, nj = int(struct["ni"]), int(struct["nj"])
    if "s_common" in struct:
        return np.repeat(np.asarray(struct["s_common"], dtype=np.float64)[None, :], ni, axis=0)
    S = np.asarray(struct["s"], dtype=np.float64)
    if S.shape != (ni, nj):
        raise StructureError(f"s の形 {S.shape} が (ni, nj) = ({ni}, {nj}) でない")
    return S


def structure_coords(struct: dict, wall_r_mm=None) -> np.ndarray:
    """構造 dict から節点座標 (N,3) [m] を作り直す。壁半径 [mm] を省くと基準値 (元の msh の壁半径)。"""
    ni, nj = int(struct["ni"]), int(struct["nj"])
    x = np.asarray(struct["x_mm"], dtype=np.float64)
    rw = np.asarray(struct["wall_r_mm"] if wall_r_mm is None else wall_r_mm, dtype=np.float64)
    R = structure_s(struct) * rw[:, None]
    out = np.zeros((ni * nj, 3))
    out[:, 0] = np.repeat(x, nj) / M_TO_MM
    out[:, 1] = R.ravel() / M_TO_MM
    return out


def write_structure_json(struct: dict, path) -> None:
    """構造 dict を JSON に書く (浮動小数は repr の 17 桁で往復一致)。"""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(struct, f, ensure_ascii=False, indent=1)
        f.write("\n")


def load_structure_json(path) -> dict:
    st = json.loads(Path(path).read_text(encoding="utf-8"))
    if st.get("format") != FORMAT or st.get("connectivity") != CONNECTIVITY:
        raise StructureError(f"構造ファイルの format/connectivity が {FORMAT}/{CONNECTIVITY} でない")
    return st


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="nozzle.msh (mesh2d の msh4.1) → mesh_structure.json")
    ap.add_argument("msh")
    ap.add_argument("--out", default=None, help="既定: msh と同じフォルダの mesh_structure.json")
    a = ap.parse_args(argv)
    st = structure_from_msh(a.msh)
    out = Path(a.out) if a.out else Path(a.msh).with_name("mesh_structure.json")
    write_structure_json(st, out)
    print(f"{out}: ni {st['ni']} nj {st['nj']} 比 {'s_common' if 's_common' in st else 's (station ごと)'} "
          f"縦線 {st['checks']['vertical_max_dev_over_L']:.2e} L")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
