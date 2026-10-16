r"""解析メッシュの構造ファイル `mesh_structure.json` (plans/accepted/tooling-nozzle-wt-deliverables.md §4.3・§6 V3)。

run の `nozzle.msh` (`meshing/mesh2d.py` の `write_msh41_2d` が書いた msh4.1) を読み、STEP の壁曲線から gmsh で
**同じ構造**のメッシュを作り直すための情報 (station $x_i$・半径方向の比 $s_{ij}$・接続の規則・物理名とタグ) を取り出す。
メッシャは再実装しない (メッシャはブランチごとに違い、今後も変わるため)。作り直しは単体スクリプト
`rebuild_from_step.py` (gmsh、forge_design に依存しない) が行う。

- 接続の規則 `mesh2d_ij_v1`: 節点番号 $(i,j)\to i\,n_j+j$ ($i$ = station、軸方向 / $j$ = 軸 0 → 壁 $n_j-1$)、
  四角形は $i$ が外側・$j$ が内側のループで $(i,j),(i{+}1,j),(i{+}1,j{+}1),(i,j{+}1)$ (反時計回り)、
  境界辺は軸 $j=0$・壁 $j=n_j-1$ ($i$ 増加の向き)、入口 $i=0$・出口 $i=n_i-1$ ($j$ 増加の向き)。
  要素番号は四角形が先で、続いて境界辺のブロックを msh に現れた順 (`write_msh41_2d` は inlet・outlet・wall・axis)。
- 縦線の検査: 各 station で $\max_j |x_{ij}-x_{i0}| \le 10^{-9} L$ ($L$ = 全長)。通れば v1 (下)。外れた格子 (壁法線の近壁層など) は
  壁の傾き $r'_i$ を渡したときだけ v2 (下) で表す。渡さなければ作らない (StructureError)。
- v1 (`forge.mesh_structure/v1`): 比 $s_{ij} = r_{ij}/r_{i,n_j-1}$。全 station の行が 1e-15 以内で同じなら `s_common` (1 行)、
  違えば `s` ($n_i\times n_j$)。msh は座標を 10 有効桁で書くので、実際の msh では行が丸めの分だけ違い、ふつうは `s` になる。
  軸の比 $s_{i0}$ は 0 に置く (軸の節点が $r=0$ であることを検査したうえで)。縦線の格子は傾きを渡しても v1 (出力は渡さないときと同じ)。
- v2 (`forge.mesh_structure/v2`、節点の規則 `wall_normal_blend_v1`; plan tooling-nozzle-card-local-mesh §5.1 #5):
  `meshing/mesh2d.py` の `wall_normal_layer` と同じ式で節点を置く。
  $(x,r)_{ij} = P_i + (1-\beta_{ij})\,d_{ij}\,(0,-1) + \beta_{ij}\,d_{ij}\,n_i$、$P_i = (x_i, r_{w,i})$、
  $n_i = (r'_i, -1)/\sqrt{1+r_i'^2}$ (内向きの壁法線)、$d_{ij} = r_{w,i}(1-s_{ij})$。
  $x_i$ は壁の節点の $x$ (軸の節点の $x$ と $10^{-12}L$ 以内で一致することを検査)。$\beta d$ と $d$ は $dX = x_{ij}-x_i$・
  $dR = r_{w,i}-r_{ij}$ から $\beta d = dX\sqrt{1+r'^2}/r'$、$d = dR + \beta d\,(1 - 1/\sqrt{1+r'^2})$ で解く
  ($|r'_i| \le 10^{-12}$ の station は $\beta = 0$)。$\beta$ は半径方向の配置 (0) と壁法線 (1) の重みなので [0, 1] に収め、
  壁の節点 ($d = 0$) は隣 ($j = n_j-2$) の値、軸の節点は 0 とする。
  v1 の項目に加えて `s` ($n_i\times n_j$、常に)・`beta` ($n_i\times n_j$)・`wall_slope` ($n_i$)・`wall_slope_source` (傾きの出典と sha256)。
  `first_layer_mm` は v1 と同じ壁の節点と $j = n_j-2$ の節点の半径の差 (壁法線の層では法線方向の厚さ × cos θ_w で厚さより小さい
  = 比の判定は厳しい側)。
  自己検査 (どちらも許容は最大 $10^{-7}$ mm かつ station の第 1 層厚の $10^{-3}$、外れたら StructureError、結果は `checks`):
  (1) 作り直し: 構造 dict ($x_i$・$r_{w,i}$・$r'_i$・$s$・$\beta$) から節点を置き直して元と比べる。$(\beta, s)$ は渡した傾きについて
      解くので、これだけでは傾きの誤りは $\beta$ が [0, 1] の外に出る節点 (ずれが $r'_i$ と逆符号の station の近壁) でしか見えない。
  (2) 規則の照合 (`rule_params` を渡したとき; 2026-10-15 主セッション判断): mesh2d の規則のパラメータ (`wall_normal_layer` の
      $[d_n, d_b]$・`wall_normal_end_stations`) から mesh2d と同じ重み
      $\beta^{rule}_{ij} = \bigl(1 - \mathrm{sm}(\mathrm{clip}((d_{ij}/r_{w,i} - d_n)/(d_b - d_n), 0, 1))\bigr)\cdot \mathrm{sm}(\mathrm{clip}(\min(i, n_i-1-i)/n_e, 0, 1))$
      ($\mathrm{sm}(t) = t^2(3-2t)$、端の断面の戻しが無ければ第 2 因子は 1) を作り、渡した傾きで節点を置き直して元と比べる。
      $\beta$ を解かずに規則から決めるので、傾きのずれは符号によらず (定数倍も) 節点のずれとして見える。パラメータは v2 の
      `rule_params` に出典と記録する。渡さなければ `checks.rule_check.status` を「未実施 (パラメータ無し)」にする (例外にしない;
      傾きの誤りの判定は STEP からの作り直しの比較 `rebuild_from_step.py --compare` だけになる)。
  **17 桁の msh が前提**: 10 有効桁の msh では近壁の $dX$ が座標の丸めに埋もれ、解いた $\beta$ が [0, 1] の外に出て (1) が、
  丸めそのもの (長さ数 m で $10^{-7}$ mm 超) で (2) が StructureError になる (case/71 の格子を 10 桁で書き直すと 9.9e-7 mm)。
  `wall_normal_layer` の格子は `mesh.msh_digits: 17` で書くこと。
- 単位: msh の座標は m、構造ファイルの長さは mm (STEP と同じ)。

usage:
  design/.venv-opt/bin/python -m forge_design.export.mesh_structure NOZZLE_MSH [--out mesh_structure.json]
  (この CLI は傾きを渡さないので縦線の格子 (v1) だけ。v2 は `deliverables.part_mesh` が run の壁ファイルから傾きを、
  prepare_info.json の格子の記録 (`mesh.params`) から規則のパラメータを作って渡す)
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

FORMAT = "forge.mesh_structure/v1"
FORMAT_V2 = "forge.mesh_structure/v2"
FORMATS = (FORMAT, FORMAT_V2)
CONNECTIVITY = "mesh2d_ij_v1"
NODE_RULE_V2 = "wall_normal_blend_v1"
VERTICAL_TOL_REL = 1e-9      # 縦線の許容 (全長に対する比; plan §4.3)
S_COMMON_TOL = 1e-15         # 比の行が全 station で同じとみなす許容
AXIS_TOL_REL = 1e-12         # 軸の節点の |r| と全節点の |z| の許容 (全長に対する比)。v2 では軸と壁の節点の x の一致にも使う
SLOPE_ZERO_TOL = 1e-12       # v2: |r′| がこれ以下の station は β = 0 (縦線; plan tooling-nozzle-card-local-mesh §5.1 #5)
SELF_CHECK_TOL_MM = 1e-7     # v2 の自己検査: 作り直した節点と元の節点の距離の最大 [mm]
SELF_CHECK_TOL_REL = 1e-3    # 同じく station の第 1 層厚に対する比 (規則の照合も同じ許容)
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
                      source=None, wall_slope=None, wall_slope_source=None, rule_params=None) -> dict:
    """msh の配列 → 構造 dict (`forge.mesh_structure/v1`、縦線でない格子は v2)。規則に従わなければ StructureError (理由つき)。

    $n_i$ = 壁辺数 + 1、$n_j$ = 入口辺数 + 1 で推定する (引数で与えたときは一致を確かめる)。
    検査: (a) 節点番号・四角形・境界辺が規則 `mesh2d_ij_v1` と完全一致、(b) 縦線、(c) 軸の節点が $r=0$
    (と全節点が $z=0$)、ほかに station が $x$ 方向に狭義増加・各 station で $r$ が $j$ について狭義増加。
    phys (物理名 → タグ) を省くと `mesh2d.PHYS` の規約を使う。
    wall_slope: station ごとの壁の傾き $r'_i = dr/dx$ (無次元、長さ $n_i$)。(b) を通る格子では使わない (v1 の出力は渡さないときと
    同じ)。(b) を通らない格子は、これがあれば v2 (節点の規則 `wall_normal_blend_v1`、モジュールの説明) を作り、無ければ StructureError。
    wall_slope_source: 傾きの出典 (dict; ファイル・sha256・求め方)。v2 の `wall_slope_source` にそのまま書く。
    rule_params: mesh2d の規則のパラメータ {"wall_normal_layer": [d_n, d_b], "wall_normal_end_stations": n_e か None,
    "source": 出典 dict (任意)}。v2 の規則の照合 (モジュールの説明) に使い、v2 の `rule_params` に書く。None なら照合は未実施。
    v1 では使わない。
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
    v2 = vmax > VERTICAL_TOL_REL
    if v2:
        if wall_slope is None:
            i = int(np.argmax(dev))
            j = int(np.argmax(np.abs(X[i] - xi[i])))
            raise StructureError(
                f"縦線の検査: station {i} (x = {xi[i] * M_TO_MM:.9g} mm) で節点 j = {j} の x が {dev[i] * M_TO_MM:.3e} mm "
                f"(= {vmax:.3e} L) ずれる (許容 {VERTICAL_TOL_REL:g} L)。格子線が縦でない (壁法線の近壁層など) 格子は v1 の構造ファイルで表せない。"
                f"壁の傾きを渡す (wall_slope: station ごとの r′ = dr/dx) と v2 (節点の規則 {NODE_RULE_V2}) で表せる")
        rp = np.asarray(wall_slope, dtype=np.float64).reshape(-1)
        if rp.shape != (ni,) or not np.all(np.isfinite(rp)):
            raise StructureError(f"壁の傾き wall_slope が長さ ni = {ni} の有限の列でない (形 {rp.shape})")
        xa = xi                                      # 軸の節点の x
        xi = X[:, -1].copy()                         # v2: station の x = 壁の節点の x
        L = float(xi[-1] - xi[0])
        if not L > 0:
            raise StructureError(f"壁の節点の全長 L = {L} ≤ 0")
        awx = float(np.abs(xa - xi).max())
        if awx > AXIS_TOL_REL * L:
            i = int(np.argmax(np.abs(xa - xi)))
            raise StructureError(f"v2 の検査: station {i} で軸の節点と壁の節点の x が {awx * M_TO_MM:.3e} mm ずれる "
                                 f"(許容 {AXIS_TOL_REL:g} L; 規則 {NODE_RULE_V2} では両端の節点が同じ x_i に乗る)")
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
    if v2:
        return _extract_v2(X, R, xi, rw, rp, L, ni, nj, bedges, phys, vmax, awx, amax, zmax, msh_sha256, source,
                           wall_slope_source, rule_params)
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


def wall_normal_blend_xr(x, rw, rp, S, B) -> tuple[np.ndarray, np.ndarray]:
    """規則 `wall_normal_blend_v1` の節点 (X, R) (ni, nj)。x・rw: station の x と壁半径 (同じ単位)、rp: 壁の傾き dr/dx、
    S・B: 比 s と重み β (ni, nj)。式は `meshing/mesh2d.py` の wall_normal_layer と同じ:
    X = x + β·d·r′/√(1 + r′²)、R = r_w − (1 − β)·d − β·d/√(1 + r′²)、d = r_w (1 − s)。"""
    x = np.asarray(x, dtype=np.float64)
    rw = np.asarray(rw, dtype=np.float64)
    rp = np.asarray(rp, dtype=np.float64)
    nrm = np.sqrt(1.0 + rp * rp)
    d = rw[:, None] * (1.0 - S)
    X = x[:, None] + B * (d * (rp / nrm)[:, None])
    R = rw[:, None] - (1.0 - B) * d - B * (d / nrm[:, None])
    return X, R


def _check_rule_params(rule_params, ni: int) -> dict:
    """規則のパラメータを検査して正規化する (mesh2d の generate_axisym_mesh と同じ条件)。不正なら StructureError。"""
    try:
        dn, db = (float(v) for v in rule_params["wall_normal_layer"])
    except (KeyError, TypeError, ValueError) as e:
        raise StructureError(f"規則のパラメータ wall_normal_layer が [d_n, d_b] でない ({rule_params!r}: {e})") from None
    if not (np.isfinite(dn) and np.isfinite(db) and 0.0 < dn < db < 1.0):
        raise StructureError(f"規則のパラメータ wall_normal_layer: 0 < d_n < d_b < 1 にする ({dn}, {db})")
    ne = rule_params.get("wall_normal_end_stations")
    if ne is not None:
        if isinstance(ne, bool) or int(ne) != ne or not (1 <= int(ne) <= (ni - 1) // 2):
            raise StructureError(f"規則のパラメータ wall_normal_end_stations: 1 以上 (ni − 1)/2 以下の整数か null にする ({ne!r}, ni {ni})")
        ne = int(ne)
    return {"wall_normal_layer": [dn, db], "wall_normal_end_stations": ne, "source": rule_params.get("source")}


def wall_normal_rule_beta(S, rule_params: dict, ni: int) -> np.ndarray:
    """mesh2d の wall_normal_layer の重み β(d) (n_i, n_j): 1 − smoothstep((d/r_w − d_n)/(d_b − d_n)) に、端の断面の戻し
    (wall_normal_end_stations = n_e) があれば smoothstep(min(i, n_i − 1 − i)/n_e) を掛ける (`meshing/mesh2d.py` と同じ式)。
    S: 比 s (d/r_w = 1 − s)。rule_params は `_check_rule_params` を通したもの。"""
    dn, db = rule_params["wall_normal_layer"]
    t = np.clip(((1.0 - S) - dn) / (db - dn), 0.0, 1.0)
    B = 1.0 - t * t * (3.0 - 2.0 * t)
    ne = rule_params["wall_normal_end_stations"]
    if ne is not None:
        ii = np.arange(ni)
        te = np.clip(np.minimum(ii, ni - 1 - ii) / float(ne), 0.0, 1.0)
        B = B * (te * te * (3.0 - 2.0 * te))[:, None]
    return B


def _solve_v2(X, R, xi, rw, rp) -> tuple[np.ndarray, np.ndarray, dict]:
    """v2: 節点 (X, R) と壁の傾き rp から比 s と重み β を解く (モジュールの説明の式)。戻り値 (S, B, 情報)。"""
    nrm = np.sqrt(1.0 + rp * rp)
    flat = np.abs(rp) <= SLOPE_ZERO_TOL
    dX = X - xi[:, None]
    dR = rw[:, None] - R
    rps = np.where(flat, 1.0, rp)
    bd = np.where(flat[:, None], 0.0, dX * (nrm / rps)[:, None])          # β·d
    d = dR + bd * (rp * rp / (nrm * (nrm + 1.0)))[:, None]                 # 1 − 1/√(1+r′²) = r′²/(√(1+r′²)(√(1+r′²)+1))
    with np.errstate(divide="ignore", invalid="ignore"):
        B = np.where(d > 0.0, bd / np.where(d > 0.0, d, 1.0), 0.0)
    B[:, -1] = B[:, -2]                                  # 壁の節点 (d = 0) は隣の値 (作り直しには効かない)
    B[:, 0] = 0.0                                        # 軸の節点は半径方向 (x = x_i、r = 0)
    out = B < 0.0
    over = B > 1.0
    info = {"beta_raw_min": float(B.min()), "beta_raw_max": float(B.max()),
            "n_beta_clamped": int(out.sum() + over.sum()), "n_flat_stations": int(flat.sum())}
    B = np.clip(B, 0.0, 1.0)                             # β は半径方向 (0) と壁法線 (1) の重み
    S = 1.0 - d / rw[:, None]
    S[:, 0] = 0.0
    S[:, -1] = 1.0
    return S, B, info


def _extract_v2(X, R, xi, rw, rp, L, ni, nj, bedges, phys, vmax, awx, amax, zmax, msh_sha256, source,
                wall_slope_source, rule_params=None) -> dict:
    """縦線でない格子の構造 dict (`forge.mesh_structure/v2`)。作った後に座標を作り直して元の節点と比べ、規則のパラメータが
    あれば規則の重みで置き直して元と比べる (自己検査 (1)(2))。"""
    rq = None if rule_params is None else _check_rule_params(rule_params, ni)
    S, B, info = _solve_v2(X, R, xi, rw, rp)
    if not np.all(np.diff(S, axis=1) > 0):
        i = int(np.argmax(~np.all(np.diff(S, axis=1) > 0, axis=1)))
        raise StructureError(f"v2: station {i} で比 s が j について狭義増加でない (壁からの距離 d が単調でない; 壁の傾きが格子と合わない可能性)")
    if phys is None:
        from forge_design.meshing.mesh2d import PHYS
        phys = dict(PHYS)
    surf = [k for k in phys if k not in SIDES]
    if not set(SIDES) <= set(phys) or len(surf) != 1:
        raise StructureError(f"物理名が {SIDES} + 面 1 つでない ({sorted(phys)})")
    physical = {nm: {"dim": 1, "tag": int(phys[nm])} for nm in SIDES}
    physical[surf[0]] = {"dim": 2, "tag": int(phys[surf[0]])}
    first = (R[:, -1] - R[:, -2]) * M_TO_MM                                # v1 と同じ半径の差 (法線方向の厚さ以下)
    st: dict = {
        "format": FORMAT_V2,
        "connectivity": CONNECTIVITY,
        "connectivity_note": "節点 (i,j) → i·nj + j (0 始まり; msh の番号は +1)。四角形 (i,j),(i+1,j),(i+1,j+1),(i,j+1) を i 外・j 内の順。"
                             "境界辺 axis j=0・wall j=nj-1 (i 増加)・inlet i=0・outlet i=ni-1 (j 増加)",
        "node_rule": NODE_RULE_V2,
        "node_rule_note": "節点 (i,j) = P_i + (1 − β_ij)·d_ij·(0, −1) + β_ij·d_ij·n_i、P_i = (x_i, r_w,i)、"
                          "n_i = (r′_i, −1)/√(1 + r′_i²) (内向きの壁法線)、d_ij = r_w,i·(1 − s_ij)。x_i = x_mm[i]、r_w,i = 壁半径 "
                          "(作り直しでは STEP の曲線の r(u_i))、r′_i = 壁の傾き dr/dx (作り直しでは STEP の曲線の (dr/du)/(dx/du))。"
                          "forge の meshing/mesh2d.py の wall_normal_layer と同じ式",
        "ni": int(ni),
        "nj": int(nj),
        "units": {"length": "mm", "msh_coords": "m", "mm_per_msh_unit": M_TO_MM,
                  "s": "1 (1 − d / r_wall、d は規則の壁からの距離)", "beta": "1 (0 = 半径方向の配置、1 = 壁法線)",
                  "wall_slope": "1 (dr/dx)"},
        "physical": physical,
        "element_tag_order": [surf[0]] + [nm for nm in bedges],
        "length_mm": L * M_TO_MM,
        "x_mm": (xi * M_TO_MM).tolist(),
        "x_note": "x_mm は station の壁の節点の x (= 軸の節点の x)。それ以外の節点の x は規則で決まる (縦線でない)",
        "s": S.tolist(),
        "beta": B.tolist(),
        "wall_slope": rp.tolist(),
        "wall_slope_source": wall_slope_source,
        "rule_params": rq,
        "rule_params_note": "mesh2d の規則のパラメータ (wall_normal_layer の [d_n, d_b]・wall_normal_end_stations) と出典。自己検査の"
                            "規則の照合に使う (作り直しは beta を使う)。null = 渡されていない (規則の照合は未実施)",
        "wall_r_mm": (rw * M_TO_MM).tolist(),
        "first_layer_mm": first.tolist(),
        "reference_note": "wall_r_mm・first_layer_mm は元の msh の各 station の壁半径と第 1 層厚 (作り直しの基準値)。第 1 層厚は v1 と同じ"
                          "壁の節点と j = nj − 2 の節点の半径の差 (壁法線の層では法線方向の厚さ × cos θ_w)",
    }
    # 自己検査: 構造 dict から作り直した節点と元の節点の距離 (plan tooling-nozzle-card-local-mesh §5.1 #5)
    Xr, Rr = wall_normal_blend_xr(np.asarray(st["x_mm"]), np.asarray(st["wall_r_mm"]), rp, S, B)
    dist = np.hypot(Xr - X * M_TO_MM, Rr - R * M_TO_MM)                     # mm
    rel = dist / first[:, None]
    k = int(np.argmax(dist))
    kr = int(np.argmax(rel))
    self_check = {"max_dist_mm": float(dist.max()), "at_node": {"i": k // nj, "j": k % nj},
                  "max_dist_over_first_layer": float(rel.max()), "ratio_at_node": {"i": kr // nj, "j": kr % nj},
                  "tol_mm": SELF_CHECK_TOL_MM, "tol_over_first_layer": SELF_CHECK_TOL_REL}
    # 規則の照合: mesh2d の規則の重み β(d) と渡した傾きで置き直す (β を解かないので傾きのずれが符号によらず見える)
    if rq is None:
        rule_check = {"status": "未実施 (パラメータ無し)",
                      "note": "規則のパラメータ (wall_normal_layer・wall_normal_end_stations) が渡されていない。傾きの誤りの判定は "
                              "STEP からの作り直しの比較 (rebuild_from_step.py --compare) だけ"}
        dq = relq = None
    else:
        Bq = wall_normal_rule_beta(S, rq, ni)
        Xq, Rq = wall_normal_blend_xr(np.asarray(st["x_mm"]), np.asarray(st["wall_r_mm"]), rp, S, Bq)
        dq = np.hypot(Xq - X * M_TO_MM, Rq - R * M_TO_MM)                   # mm
        relq = dq / first[:, None]
        kq, kqr = int(np.argmax(dq)), int(np.argmax(relq))
        # β の差 (情報): 壁の節点 (d = 0) と傾きの無い station (|r′| ≤ slope_zero_tol、解いた β は 0) は β が節点に効かないので除く
        dB = np.where((np.abs(rp) <= SLOPE_ZERO_TOL)[:, None], 0.0, np.abs(Bq - B))[:, :-1]
        kb = int(np.argmax(dB))
        rule_check = {"status": "実施", "max_dist_mm": float(dq.max()), "at_node": {"i": kq // nj, "j": kq % nj},
                      "max_dist_over_first_layer": float(relq.max()), "ratio_at_node": {"i": kqr // nj, "j": kqr % nj},
                      "beta_rule_minus_solved_max_abs": float(dB.max()),
                      "beta_at_node": {"i": kb // (nj - 1), "j": kb % (nj - 1)},
                      "tol_mm": SELF_CHECK_TOL_MM, "tol_over_first_layer": SELF_CHECK_TOL_REL,
                      "note": "mesh2d の規則の重み β(d) (rule_params) と wall_slope で置き直した節点と元の節点の距離。"
                              "beta_rule_minus_solved は情報 (壁の節点と傾きの無い station (|r′| ≤ slope_zero_tol) を除く)"}
    st["checks"] = {"vertical_max_dev_over_L": vmax, "vertical_tol_over_L": VERTICAL_TOL_REL,
                    "axis_wall_x_max_dev_over_L": awx / L, "axis_max_abs_r_m": amax, "max_abs_z_m": zmax,
                    "slope_zero_tol": SLOPE_ZERO_TOL, **info, "self_check": self_check,
                    "self_check_note": "(β, s) は渡した傾きについて解くので、作り直しの一致だけでは傾きの正しさを示さない (傾きの誤りが見えるのは "
                                       "β が [0, 1] の外に出る節点だけ = n_beta_clamped)。傾きは rule_check (規則の照合) と STEP からの"
                                       "作り直しの比較で確かめる",
                    "rule_check": rule_check}
    st["source"] = {"msh": None if source is None else str(source), "sha256": msh_sha256}
    if not (dist.max() <= SELF_CHECK_TOL_MM and rel.max() <= SELF_CHECK_TOL_REL):
        raise StructureError(
            f"v2 の自己検査: 作り直した節点が元から最大 {dist.max():.3e} mm (station {k // nj}・j = {k % nj})、第 1 層厚の "
            f"{rel.max():.3e} 倍 (station {kr // nj}・j = {kr % nj}) ずれる (許容 {SELF_CHECK_TOL_MM:g} mm かつ {SELF_CHECK_TOL_REL:g}; "
            f"β が [0, 1] の外に出た節点 {info['n_beta_clamped']}、解いた β の範囲 [{info['beta_raw_min']:.9g}, "
            f"{info['beta_raw_max']:.9g}])。壁の傾き wall_slope が格子と合わない可能性")
    if dq is not None and not (dq.max() <= SELF_CHECK_TOL_MM and relq.max() <= SELF_CHECK_TOL_REL):
        rc = rule_check
        raise StructureError(
            f"v2 の規則の照合: mesh2d の規則 (wall_normal_layer {rq['wall_normal_layer']}・wall_normal_end_stations "
            f"{rq['wall_normal_end_stations']}) の重みと渡した傾きで置き直した節点が元から最大 {rc['max_dist_mm']:.3e} mm "
            f"(station {rc['at_node']['i']}・j = {rc['at_node']['j']})、第 1 層厚の {rc['max_dist_over_first_layer']:.3e} 倍 "
            f"(station {rc['ratio_at_node']['i']}・j = {rc['ratio_at_node']['j']}) ずれる (許容 {SELF_CHECK_TOL_MM:g} mm かつ "
            f"{SELF_CHECK_TOL_REL:g}; |β_規則 − β_解| 最大 {rc['beta_rule_minus_solved_max_abs']:.3e})。壁の傾き wall_slope か"
            "規則のパラメータが格子と合わない")
    return st


def structure_beta(struct: dict) -> np.ndarray:
    """v2 の構造 dict の重み β (n_i, n_j)。"""
    ni, nj = int(struct["ni"]), int(struct["nj"])
    B = np.asarray(struct["beta"], dtype=np.float64)
    if B.shape != (ni, nj):
        raise StructureError(f"beta の形 {B.shape} が (ni, nj) = ({ni}, {nj}) でない")
    return B


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


def structure_coords(struct: dict, wall_r_mm=None, wall_slope=None) -> np.ndarray:
    """構造 dict から節点座標 (N,3) [m] を作り直す。壁半径 [mm] を省くと基準値 (元の msh の壁半径)。
    v2 は規則 `wall_normal_blend_v1` で置く (壁の傾きを省くと構造 dict の `wall_slope`)。"""
    ni, nj = int(struct["ni"]), int(struct["nj"])
    x = np.asarray(struct["x_mm"], dtype=np.float64)
    rw = np.asarray(struct["wall_r_mm"] if wall_r_mm is None else wall_r_mm, dtype=np.float64)
    if struct.get("format") == FORMAT_V2:
        rp = np.asarray(struct["wall_slope"] if wall_slope is None else wall_slope, dtype=np.float64)
        X, R = wall_normal_blend_xr(x, rw, rp, structure_s(struct), structure_beta(struct))
        out = np.zeros((ni * nj, 3))
        out[:, 0] = X.ravel() / M_TO_MM
        out[:, 1] = R.ravel() / M_TO_MM
        return out
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
    if st.get("format") not in FORMATS or st.get("connectivity") != CONNECTIVITY:
        raise StructureError(f"構造ファイルの format/connectivity が {' か '.join(FORMATS)}/{CONNECTIVITY} でない")
    if st["format"] == FORMAT_V2 and st.get("node_rule") != NODE_RULE_V2:
        raise StructureError(f"v2 の構造ファイルの節点の規則 node_rule が {NODE_RULE_V2} でない ({st.get('node_rule')!r})")
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
