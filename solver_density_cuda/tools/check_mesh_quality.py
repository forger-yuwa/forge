#!/usr/bin/env python3
"""forge 入力 H5 メッシュの品質指標 (アスペクト比・スキューネス) を計算・判定する。

使い方:
    python3 check_mesh_quality.py mesh.h5 [--ar-max 1000] [--skew-max 0.9] [--mode auto|2d|3d]

- アスペクト比 (AR): 各セルの 最長辺/最短辺。薄い境界層セルの監視用。
- スキューネス (equiangle skew): 内角 θ から max[(θmax-ideal)/(180-ideal),(ideal-θmin)/ideal]
  (0=直交/正三角、1=退化)。四角形 ideal=90、三角形 ideal=60。
- 既定しきい値: AR ≤ 1000, skew ≤ 0.9 (本リポジトリの計算前メッシュ品質ルール、
  AGENTS.md / forge-calculation-workflow.md)。VERDICT を PASS/FAIL で返す。

2D / 3D の扱い (重要):
- **2D モード** (準2D・z 押し出し・軸対称): セルを x-y 平面に射影し多角形として測る (従来挙動)。
- **3D モード** (完全非構造 3D: tet/prism/hex/pyramid が z 方向にも構造を持つ): x-y 射影は
  無意味 (z 平行エッジが一点に潰れ AR=inf 等の偽値) なため、**各セルの面ごとに equiangle skew
  を測りその最大**を skew、**全辺の最長/最短**を AR とする (真の 3D 形状で判定)。
- `--mode auto` (既定): 3D セル型を含み z 平面が 3 枚以上なら 3D、それ以外は 2D と自動判定。

`solver_density_cuda/tools/res_h5_to_vtu.py` の parse_conne を流用して CONNE を読む。
node (median-dual) 変換の h5 は `/MESH/CONNE` が双対 CV (多角形/多面体) なので、`/VIZMESH/CONNE` (primal
セル) があればそちらを測る (品質は primal の性質。cell 変換をやり直す必要はない)。
計算は VTK セル型ごとに numpy でベクトル化 (2026-09-08: 2M セル hex で 15 分 → 数秒)。`--loop` で旧の
セル毎ループ (検証用) を使う。

壁直交性と AR の緩和条件 (`--bcond bcondConfig.yaml`、plan tooling-generator-defaults-and-launch-gate §4.4・§5.1 #8・#8c、2026-10-14):
- 壁の bcond (`wall`・`wall_isothermal`) の各境界面 (2D は辺) について、隣接セルの「壁上の節点から壁の外の節点へ向かう辺」
  (横断辺) と壁法線のなす角 φ を測る。四角形・六面体・プリズム (三角形の面から入る層) が対象で、三角形・四面体・ピラミッドと
  四角形の面から入るプリズムは「対象外」として数だけ出す。境界面は node 変換なら `/BCONDS/<physID>/vizBfaceNodes`、
  cell 変換なら `iPlanes` → `/PLANES/STRUCT` の primal 節点で、bcond の physID と対応づける。
- 壁から層を数える (列): 壁の面 1 つを列の始まりとし、入ったセルの向かいの面を渡って隣へ進む (層 1, 2, ...)。層のセルの φ は
  「入った面の法線と横断辺のなす角」、内角は横断辺を含む側面の角の 90° からのずれ。
- **AR の緩和は列ごと** (#8c): 1000 < AR <= `--ar-max` (最大 5000) のセルは、壁の面からそのセルまでの同じ列の全層
  (層 1..そのセル) で φ <= `--phi-construct` かつ内角のずれ <= `--theta-tol` のときだけ認める。セルは最も近い (層の番号が最小の)
  列で判定し、同じ層の番号で 2 本以上の列が届くときは全部の列で満たすこと。確かめられないセル (どの列も届かない・層でない型・
  退化) は認めない。**physID は報告の単位で、合否の単位にしない** (同じ座標・接続・壁種なら、壁の physID の分け方で緩和の可否・
  警告の対象は変わらない)。
- **`--bcond` のとき AR の外れ値に SOFT-PASS を当てない** (#8c): AR > `--ar-max` のセル (AR > 5000 を含む) と、緩和の対象で
  確かめられないセルは、件数によらず FAIL。skew の外れ値 (0.1 % 未満) の SOFT-PASS は従来どおり。`--ar-max` > 5000 はエラー。
- **閾値は構成の判定であって精度の限界ではない**: `--phi-construct`・`--theta-tol` は「層が壁法線に沿って作られたか」を
  分けるための暫定値 (半径方向の線の 40° 傾斜壁 φ = 40° と壁法線の層 φ ≈ 0 の間)。どこまで傾いてよいかの精度の予算は
  壁応力・熱流束・δ* の A/B で別に決める (plan §4.4・§6.4)。
- AR <= 1000 の近壁 (層 1) で φ > `--phi-construct` のセルは WARN (数・最大 φ・位置)。FAIL にしない。角部のセルも含める。
- 角部 = 壁の面の (領域の内側へ向けた) 法線が `--corner-angle` 以上折れる節点 (幾何の角。45° は分類の暫定値)。境界のラベルの
  継ぎ目 (壁の physID が変わるだけで折れが小さい節点) は角部にしない。角部は報告 (節点の数・最大の折れ角と位置・角部のセルの
  φ の最大と位置) を分けるだけで、緩和の判定・近壁の WARN からは外さない (#8c)。
- 結果の再利用: 鍵 (メッシュの不変部分のハッシュ・physID と壁の面の集合の対応・閾値・検査器の版) が一致したときだけ
  `<h5>.quality.json` の結果を使う。再利用したかどうかを出力に書く (`--no-cache` で使わない)。
- `--bcond` を渡さないときの VERDICT・終了コードは従来どおり (壁直交性は測らない)。その判定は投入の関門の合格の根拠に
  ならない旨を必ず 1 行出す。`--ar-max` を 1000 より上げると緩和の条件を確かめずに通す旧判定なので、その旨も注意として出す。
"""
import sys, os, argparse, hashlib, json
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne

# VTK cell type -> 面 (セル局所節点 index)。parse_conne の節点順 (gmsh/XDMF/VTK) 準拠。
# TETRA=10, HEXAHEDRON=12, WEDGE(prism)=13, PYRAMID=14
FACES_3D = {
    10: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
    12: [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    13: [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    14: [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
VTK_3D = set(FACES_3D)


def poly_skew(pts):
    """多角形 (n,2 or n,3) の equiangle skew を返す (cell_metrics の skew 部と同一定義)。"""
    n = len(pts)
    angs = []
    for i in range(n):
        a = pts[(i - 1) % n] - pts[i]
        b = pts[(i + 1) % n] - pts[i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na == 0 or nb == 0:
            angs.append(0.0); continue
        c = np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0)
        angs.append(np.degrees(np.arccos(c)))
    angs = np.array(angs)
    ideal = 90.0 if n == 4 else (60.0 if n == 3 else 180.0 * (n - 2) / n)
    return max((angs.max() - ideal) / (180.0 - ideal), (ideal - angs.min()) / ideal)


def cell_metrics(pts):
    """pts: (n,2) 多角形頂点 (CCW or CW)。(aspect_ratio, skew) を返す (2D モード用)。"""
    n = len(pts)
    edges = [np.linalg.norm(pts[(i + 1) % n] - pts[i]) for i in range(n)]
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    return ar, poly_skew(pts)


# VTK セル型ごとの四面体分割 (符号付き体積の合計を取るため)
TETRA_OF = {
    10: [(0, 1, 2, 3)],                                              # tetra
    12: [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (1, 4, 5, 6), (3, 4, 6, 7)],  # hexa
    13: [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)],                  # wedge (prism)
    14: [(0, 1, 2, 4), (0, 2, 3, 4)],                                # pyramid
}


# 六面体の各頂点で、その頂点から出る 3 辺 (右手系の順) が張る体積 = 頂点 Jacobian。VTK/gmsh の頂点順
HEX_CORNER = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7), 4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}


def hex_corner_jacobians(coord3, hx):
    """(n, 8) の頂点 Jacobian。全部同符号なら妥当な六面体、符号が混じれば局所反転。"""
    P = np.asarray(coord3, dtype=np.float64)
    return np.stack([np.einsum("ij,ij->i", np.cross(P[hx[:, a]] - P[hx[:, c]], P[hx[:, b]] - P[hx[:, c]]), P[hx[:, d]] - P[hx[:, c]])
                     for c, (a, b, d) in HEX_CORNER.items()], 1)


def cell_volumes(coord3, conn, offs, types):
    """セルの符号付き体積。型が未知なら NaN。**成立性検査専用**。戻り値 (vols, 局所反転した六面体の数)。

    六面体は四面体分割の和を使わない (2026-09-22): 5 四面体分割は向かい合う面を別の対角線で切るので、
    **面の反り (曲面壁の弦の矢高) が層厚より大きい薄い壁層セル**では分割どうしが交差して和が負になり、頂点 Jacobian が
    8 つとも正の妥当なセルを「向きの不整合」と誤判定した (case/46 接続模型のフィレット壁層 111 セル)。逆に和だけでは
    局所反転 (最小頂点 Jacobian < 0) を見逃す (codex 2026-09-21)。六面体は**頂点 Jacobian の平均を体積、符号の混在を局所反転**とする。
    """
    vols = np.full(len(offs), np.nan)
    offs = np.asarray(offs); types = np.asarray(types); conn = np.asarray(conn)
    starts = np.r_[0, offs[:-1]]; ishex = (types == 12) & ((offs - starts) == 8)
    n_inverted = 0
    if ishex.any():
        hx = conn[(starts[ishex][:, None] + np.arange(8)[None, :])]
        for i0 in range(0, len(hx), 500000):
            J = hex_corner_jacobians(coord3, hx[i0:i0 + 500000])
            v = J.mean(1); mixed = (J.min(1) * J.max(1) <= 0.0) & np.isfinite(v)
            n_inverted += int(mixed.sum())
            vols[np.flatnonzero(ishex)[i0:i0 + 500000]] = v
    s = 0
    for i, o in enumerate(offs):
        idx = conn[s:o]; s = o
        if ishex[i]:
            continue
        tets = TETRA_OF.get(int(types[i]))
        if tets is None or len(idx) < 4:
            continue
        pts = coord3[idx]
        v = 0.0
        okc = True
        for t in tets:
            if max(t) >= len(pts):
                okc = False
                break
            p0, p1, p2, p3 = pts[t[0]], pts[t[1]], pts[t[2]], pts[t[3]]
            v += np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0
        vols[i] = v if okc else np.nan
    return vols, n_inverted


def cell_metrics_3d(pts, vtk):
    """pts: (m,3) セル節点。面ごとの equiangle skew の最大と、全辺の 最長/最短 を返す。"""
    faces = FACES_3D[vtk]
    skew = 0.0
    edges = []
    for fc in faces:
        fp = pts[list(fc)]
        skew = max(skew, poly_skew(fp))
        for i in range(len(fc)):
            edges.append(np.linalg.norm(fp[(i + 1) % len(fc)] - fp[i]))
    emin, emax = min(edges), max(edges)
    ar = emax / emin if emin > 0 else float("inf")
    return ar, skew


def _poly_skew_batch(P):
    """P: (nc, n, d) 多角形頂点列。各多角形の equiangle skew (nc,) をベクトル計算 (poly_skew と同一定義)。"""
    n = P.shape[1]
    a = np.roll(P, 1, axis=1) - P          # pts[i-1] - pts[i]
    b = np.roll(P, -1, axis=1) - P         # pts[i+1] - pts[i]
    na = np.linalg.norm(a, axis=2); nb = np.linalg.norm(b, axis=2)
    den = na * nb
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(den > 0, np.einsum("ijk,ijk->ij", a, b) / np.where(den > 0, den, 1.0), 1.0)
    ang = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    ang = np.where(den > 0, ang, 0.0)      # 退化辺 (長さ 0) は角 0 (poly_skew と同じ)
    ideal = 90.0 if n == 4 else (60.0 if n == 3 else 180.0 * (n - 2) / n)
    return np.maximum((ang.max(axis=1) - ideal) / (180.0 - ideal), (ideal - ang.min(axis=1)) / ideal)


def _edge_ar_batch(P):
    """P: (nc, n, d)。各多角形の辺長 (nc, n) を返す。"""
    return np.linalg.norm(np.roll(P, -1, axis=1) - P, axis=2)


def metrics_vectorized(coord, conn, offs, types, is3d):
    """セル型ごとにまとめて (ars, skews) を計算する (cell_metrics / cell_metrics_3d と同一定義)。"""
    ncells = len(offs)
    ars = np.empty(ncells); skews = np.empty(ncells)
    starts = np.concatenate(([0], offs[:-1]))
    sizes = offs - starts
    for t in np.unique(types):
        sel = np.where(types == t)[0]
        nn = np.unique(sizes[sel])
        for m in nn:                                   # 同型でも節点数が違えば分ける (安全側)
            idx = sel[sizes[sel] == m]
            gather = starts[idx][:, None] + np.arange(m)[None, :]
            P = coord[conn[gather]]                     # (nc, m, d)
            if not is3d:
                e = _edge_ar_batch(P)
                emin = e.min(axis=1); emax = e.max(axis=1)
                ars[idx] = np.where(emin > 0, emax / np.where(emin > 0, emin, 1.0), np.inf)
                skews[idx] = _poly_skew_batch(P)
            else:
                faces = FACES_3D.get(int(t))
                if faces is None:
                    raise SystemExit(f"unsupported 3D VTK type {t}")
                sk = np.zeros(len(idx)); emin = np.full(len(idx), np.inf); emax = np.zeros(len(idx))
                for fc in faces:
                    fp = P[:, list(fc), :]
                    sk = np.maximum(sk, _poly_skew_batch(fp))
                    e = _edge_ar_batch(fp)
                    emin = np.minimum(emin, e.min(axis=1)); emax = np.maximum(emax, e.max(axis=1))
                ars[idx] = np.where(emin > 0, emax / np.where(emin > 0, emin, 1.0), np.inf)
                skews[idx] = sk
    return ars, skews


# ---------------------------------------------------------------------------------------------------------------------------
# 壁直交性と AR の緩和条件 (plan tooling-generator-defaults-and-launch-gate §4.4・§5.1 #8、2026-10-14)
# ---------------------------------------------------------------------------------------------------------------------------
CHECKER_VERSION = "wall-ortho-v2"   # v2 (#8c): 列ごとの全層・幾何の角だけを角部・--bcond で AR の外れ値に SOFT-PASS なし
AR_BASE = 1000.0               # AGENTS.md の既定の上限。これを超え --ar-max 以下のセルが「緩和の対象」
AR_RELAX_CAP = 5000.0          # 緩和の上限 (AGENTS.md の例外、2026-09-12 ユーザ決定)
WALL_KINDS = ("wall", "wall_isothermal")
# 構成の判定の暫定値 (精度の限界ではない; plan §4.4・§5.1 #8 の記録)。2026-10-14 の測定 (格納 float64):
#   - 半径方向の線の 40° 傾斜壁 (合成): φ = 40.00°。半径方向の線は傾斜 θ_w の壁で φ = θ_w。
#   - case/45 の生産の格子 (problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml、wall_normal_layer [0.02, 0.8]、
#     ni 4719・nj 160): 層 1 の φ 最大 0.06°、AR > 1000 の層 1..19 で φ 最大 0.07°・内角のずれ最大 0.07°。
#     同じ壁で wall_normal_layer なし (半径方向の線) にすると φ 最大 26.6° (壁の傾き 26.6°)。
#   5° は生産の 0.07° の約 70 倍、傾き 5° を超える壁の半径方向の線は拒む。曲面壁で節点ごとの法線に置いた線は
#   φ ≈ (1 セルあたりの壁の折れ角)/2 なので、折れ角 10°/セルまでは通る。内角のずれも同じ値 (生産で 0.07°)。
PHI_CONSTRUCT_DEFAULT = 5.0    # [deg]
THETA_TOL_DEFAULT = 5.0        # [deg] 緩和の対象の層の内角の 90° からのずれ
CORNER_ANGLE_DEFAULT = 45.0    # [deg] 壁の面の法線がこれ以上折れる節点を角部とする (報告を分けるための分類の暫定値)
PHI_BINS = (1.0, 5.0, 10.0, 20.0)

# 2D (平面の三角形・四角形) の局所の辺。3D は FACES_3D。
FACES_2D = {5: [(0, 1), (1, 2), (2, 0)], 9: [(0, 1), (1, 2), (2, 3), (3, 0)]}
# 層として進める (入った面 → 向かいの面)。三角形・四面体・ピラミッドと、四角形の面から入るプリズムは層でない
OPPOSITE = {9: {0: 2, 1: 3, 2: 0, 3: 1}, 12: {0: 1, 1: 0, 2: 4, 4: 2, 3: 5, 5: 3}, 13: {0: 1, 1: 0}}
CELL_EDGES = {
    9: [(0, 1), (1, 2), (2, 3), (3, 0)],
    12: [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
    13: [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
}
VTK_NAME = {5: "三角形", 9: "四角形", 10: "四面体", 12: "六面体", 13: "プリズム", 14: "ピラミッド"}


def _local_faces(vtk, topo3d):
    return FACES_3D.get(vtk, []) if topo3d else FACES_2D.get(vtk, [])


def _layer_geometry(vtk, lf, topo3d):
    """(型, 入った面) → (入った面の局所節点, 横断辺 [(壁側, 外側)], 側面の多角形)。層として進めないなら None。"""
    if lf not in OPPOSITE.get(vtk, {}) or (vtk == 9) == topo3d:
        return None
    faces = _local_faces(vtk, topo3d)
    ent, opp = faces[lf], faces[OPPOSITE[vtk][lf]]
    trans = []
    for u, v in CELL_EDGES[vtk]:
        if u in ent and v in opp:
            trans.append((u, v))
        elif v in ent and u in opp:
            trans.append((v, u))
    if vtk == 9:
        lat = [(0, 1, 2, 3)]                 # 2D の四角形は自身が側面 (横断辺と壁に沿う辺の角)
    else:
        lat = [f for i, f in enumerate(faces) if i not in (lf, OPPOSITE[vtk][lf])]
    return ent, trans, lat


def _face_normals(P):
    """P: (n, k, 3)。k == 2 は x-y 面内の辺の法線 (2D)、k >= 3 は Newell の法線。正規化しない。"""
    if P.shape[1] == 2:
        t = P[:, 1] - P[:, 0]
        return np.stack([t[:, 1], -t[:, 0], np.zeros(len(t))], 1)
    return np.cross(P, np.roll(P, -1, axis=1)).sum(axis=1)


def _phi_deg(E, n):
    """E: (n, m, 3) 辺、n: (n, 3) 面の法線 → (n, m) のなす角 [deg] (0 = 法線に沿う)。退化は NaN。"""
    ne = np.linalg.norm(E, axis=2)
    nn = np.linalg.norm(n, axis=1)[:, None]
    den = ne * nn
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.abs(np.einsum("ijk,ik->ij", E, n)) / den
    return np.where(den > 0, np.degrees(np.arccos(np.clip(c, 0.0, 1.0))), np.nan)


def _poly_angles(P):
    """P: (n, k, 3) 多角形 → (n, k) の内角 [deg]。退化は NaN。"""
    a = np.roll(P, 1, axis=1) - P
    b = np.roll(P, -1, axis=1) - P
    den = np.linalg.norm(a, axis=2) * np.linalg.norm(b, axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.einsum("ijk,ijk->ij", a, b) / den
    return np.where(den > 0, np.degrees(np.arccos(np.clip(c, -1.0, 1.0))), np.nan)


def _round_allow_deg(Ea_len, ra, rb):
    """辺 (長さ Ea_len、両端の丸めの上限 ra, rb) の向きの丸めの上限 [deg] (小角近似、90° で打ち切り)。"""
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where(Ea_len > 0, (ra + rb) / Ea_len, np.inf)
    return np.degrees(np.minimum(v, np.pi / 2))


def _face_keys(nodes, nnodes):
    """(n, k <= 4) の節点 → 向きと順序によらない面の鍵 (k1, k2)。k1 は大きい 2 節点 (2D の辺でも一意になる側)、
    k2 は残り (パディング -1 を含む)。主キー k1 で二分探索するので、k1 が区別できない並べ方にしない。"""
    s = np.full((len(nodes), 4), -1, dtype=np.int64)
    s[:, :nodes.shape[1]] = nodes
    s.sort(axis=1)
    m = np.int64(nnodes + 1)
    return (s[:, 2] + 1) * m + (s[:, 3] + 1), (s[:, 0] + 1) * m + (s[:, 1] + 1)


def read_bcond_kinds(path):
    """bcondConfig.yaml → {physID: (名前, kind)}。読めなければ例外。"""
    import yaml
    with open(path) as fh:
        y = yaml.safe_load(fh) or {}
    out = {}
    for nm, d in y.items():
        if isinstance(d, dict) and "physID" in d and "kind" in d:
            out[int(d["physID"])] = (str(nm), str(d["kind"]))
    if not out:
        raise ValueError("physID と kind を持つ境界が 1 つも無い")
    return out


def _planes_offsets(st):
    """cell 変換の /PLANES/STRUCT ([nNodes, nodes..., nCells, cells...] の繰り返し) → 各面の先頭 (nNodes の位置)。
    全体を Python のリストにしない (3D の数百万面で GB 単位になる)。"""
    offs, i, n = [], 0, len(st)
    while i < n:
        offs.append(i)
        i += 1 + int(st[i])
        i += 1 + int(st[i])
    return np.asarray(offs, dtype=np.int64)


def read_boundary_faces(f, pids):
    """{physID: [節点の配列, ...]} (primal の節点)。node 変換は vizBfaceNodes、cell 変換は iPlanes → PLANES/STRUCT。
    取れない physID は None。"""
    out, planes = {}, None
    for pid in pids:
        key = "BCONDS/%d" % pid
        if key not in f:
            out[pid] = None
            continue
        g = f[key]
        if "vizBfaceNodes" in g and "vizBfaceSizes" in g:
            nodes = np.asarray(g["vizBfaceNodes"], dtype=np.int64)
            sizes = np.asarray(g["vizBfaceSizes"], dtype=np.int64)
            out[pid] = np.split(nodes, np.cumsum(sizes)[:-1]) if len(sizes) else []
        elif "iPlanes" in g and "PLANES/STRUCT" in f and "VIZMESH" not in f:
            if planes is None:
                st = np.asarray(f["PLANES/STRUCT"], dtype=np.int64)
                planes = (st, _planes_offsets(st))
            st, po = planes
            out[pid] = [st[po[i] + 1:po[i] + 1 + st[po[i]]] for i in np.asarray(g["iPlanes"], dtype=np.int64)]
        else:
            out[pid] = None
    return out


def _sha_update_ds(h, f, name):
    if name in f:
        d = np.ascontiguousarray(np.asarray(f[name]))
        h.update(f"{name}|{d.dtype.str}|{d.shape}|".encode())
        h.update(d.tobytes())


def mesh_invariant_sha(f):
    """メッシュの不変部分 (座標・接続・境界の所属) のハッシュ。/VALUE (restart で変わる) は含めない。"""
    h = hashlib.sha256()
    for name in ("MESH/COORD", "MESH/CONNE", "VIZMESH/CONNE"):
        _sha_update_ds(h, f, name)
    if "BCONDS" in f:
        for k in sorted(f["BCONDS"], key=lambda x: int(x) if str(x).isdigit() else 0):
            for ds in ("iPlanes", "iBPlanes", "iCells", "vizBfaceNodes", "vizBfaceSizes"):
                _sha_update_ds(h, f, "BCONDS/%s/%s" % (k, ds))
    return h.hexdigest()


def _faces_sha(faces):
    """壁の面の集合 (向き・順序によらない) のハッシュ。"""
    keys = sorted(tuple(sorted(int(x) for x in fc)) for fc in faces)
    return hashlib.sha256(json.dumps(keys).encode()).hexdigest()


def checker_sha():
    h = hashlib.sha256()
    for fn in (os.path.abspath(__file__), os.path.join(os.path.dirname(os.path.abspath(__file__)), "res_h5_to_vtu.py")):
        with open(fn, "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()


def node_rounding_radius(coord3):
    """節点の位置の丸め誤差の上限 [長さ] (格納の dtype の ulp の半分 × √3)。float32 の h5 では薄い層の辺の向きが
    ulp / 層厚 [rad] だけ揺れる (AR 4000・座標 0.5 m・層厚 1 µm で約 6°)。構成の判定はこれを差し引いて行う。"""
    c = np.abs(np.asarray(coord3))
    ulp = np.spacing(c.astype(c.dtype if np.issubdtype(c.dtype, np.floating) else np.float64))
    return (np.sqrt(3.0) * 0.5 * ulp.max(axis=1)).astype(np.float64)


def wall_orthogonality(coord3, conn, offs, types, walls, ars, a, rnode=None):
    """壁直交性・層・AR の緩和の判定。walls: {physID: (名前, kind, [面の節点の配列])} (壁だけ)。
    rnode: 節点の丸め誤差の上限 (node_rounding_radius)。**判定は格納された座標そのもの (ソルバが計算に使う形) の φ・内角で行う**。
    丸めの上限は、条件を満たさない理由が格納精度の揺れの範囲か (丸めを差し引けば満たすか) を分けて出すためだけに使う。

    AR の緩和と近壁の WARN は列ごと (壁の面 1 つから向かいの面を渡る層の並び) に判定し、physID は報告の単位にだけ使う
    (#8c: 同じ座標・接続・壁種なら、壁の physID の分け方で緩和の可否・警告の対象が変わらない)。

    戻り dict: lines (出力行)・n_unverified (緩和の対象で確かめられないセル)・warn (近壁の φ の WARN 数)・
    indeterminate (判定不能の理由のリスト)・summary (JSON 用)・masks (セルごとの対象・確認・WARN・角部・層の番号 (最も近い列の層。行進を打ち切った先は 0)。試験用)。"""
    lines, indet = [], []
    say = lines.append
    conn = np.asarray(conn, dtype=np.int64)
    offs = np.asarray(offs, dtype=np.int64)
    types = np.asarray(types).astype(np.int64)
    ncells, nnodes = len(offs), len(coord3)
    starts = np.r_[0, offs[:-1]]
    P3 = np.asarray(coord3, dtype=np.float64)
    R3 = np.zeros(nnodes) if rnode is None else np.asarray(rnode, dtype=np.float64)
    topo3d = bool(np.isin(types, list(VTK_3D)).any())
    maxlf = 6 if topo3d else 4

    # --- 面の表 (セルの局所の面 → 隣のセル・隣での局所の面) ---
    K1, K2, FC, FL = [], [], [], []
    for t in np.unique(types):
        faces = _local_faces(int(t), topo3d)
        sel = np.flatnonzero(types == t)
        for lf, fc in enumerate(faces):
            k1, k2 = _face_keys(conn[starts[sel][:, None] + np.asarray(fc)[None, :]], nnodes)
            K1.append(k1); K2.append(k2); FC.append(sel); FL.append(np.full(len(sel), lf, dtype=np.int64))
    K1 = np.concatenate(K1); K2 = np.concatenate(K2); FC = np.concatenate(FC); FL = np.concatenate(FL)
    order = np.lexsort((K2, K1))
    s1, s2 = K1[order], K2[order]
    nbr_c = np.full((ncells, maxlf), -1, dtype=np.int64)
    nbr_l = np.full((ncells, maxlf), -1, dtype=np.int64)
    same = np.flatnonzero((s1[1:] == s1[:-1]) & (s2[1:] == s2[:-1]))
    ia, ib = order[same], order[same + 1]
    nbr_c[FC[ia], FL[ia]] = FC[ib]; nbr_l[FC[ia], FL[ia]] = FL[ib]
    nbr_c[FC[ib], FL[ib]] = FC[ia]; nbr_l[FC[ib], FL[ib]] = FL[ia]
    opp_tab = np.full((16, maxlf), -1, dtype=np.int64)
    for t, mp in OPPOSITE.items():
        if (t == 9) != topo3d:
            for i, j in mp.items():
                opp_tab[t, i] = j

    def lookup(nodes_list):
        """境界面 → (セル, 局所の面)。見つからない・2 つ以上のセルに属する面は -1。"""
        if not nodes_list:
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        kk = max(len(x) for x in nodes_list)
        arr = np.full((len(nodes_list), kk), -1, dtype=np.int64)
        for i, x in enumerate(nodes_list):
            arr[i, :len(x)] = x
        q1, q2 = _face_keys(arr, nnodes)
        # パディング (-1) は _face_keys の鍵で 0 番の扱いになり、局所の面の鍵 (-1 パディング) と一致する
        lo = np.searchsorted(s1, q1, "left"); hi = np.searchsorted(s1, q1, "right")
        cell = np.full(len(q1), -1, dtype=np.int64); lfo = np.full(len(q1), -1, dtype=np.int64)
        for i in range(len(q1)):
            hit = [order[j] for j in range(lo[i], hi[i]) if s2[j] == q2[i]]
            if len(hit) == 1:
                cell[i], lfo[i] = FC[hit[0]], FL[hit[0]]
        return cell, lfo

    # --- 層ごとの φ・内角 (グループ = (型, 入った面)) ---
    def layer_metrics(cells, lfs):
        """戻り (φ, 内角のずれ, φ の下限, 内角のずれの下限, 層の型か)。下限は丸めの上限を差し引いた値。"""
        phi = np.full(len(cells), np.nan); ang = np.full(len(cells), np.nan)
        phl = np.full(len(cells), np.nan); anl = np.full(len(cells), np.nan)
        is_layer = np.zeros(len(cells), dtype=bool)
        if not len(cells):
            return phi, ang, phl, anl, is_layer
        code = types[cells] * 16 + lfs
        for cval in np.unique(code):
            t, lf = int(cval // 16), int(cval % 16)
            g = _layer_geometry(t, lf, topo3d)
            if g is None:
                continue
            ent, trans, lat = g
            idx = np.flatnonzero(code == cval)
            nn = offs[cells[idx]] - starts[cells[idx]]
            ids = conn[starts[cells[idx]][:, None] + np.arange(int(nn[0]))[None, :]]
            P, R = P3[ids], R3[ids]

            def rowmax(x):
                return np.where(np.isnan(x).any(axis=1), np.nan, np.nan_to_num(x, nan=0.0).max(axis=1))

            ph, pl = _phi_parts(P, R, ent, trans)
            phi[idx], phl[idx] = rowmax(ph), rowmax(pl)
            dev, dlo = np.zeros(len(idx)), np.zeros(len(idx))
            for poly in lat:
                dv, dl = _angle_dev_parts(P[:, list(poly)], R[:, list(poly)])
                dev = np.where(np.isnan(dev), np.nan, np.maximum(dev, rowmax(dv)))
                dlo = np.where(np.isnan(dlo), np.nan, np.maximum(dlo, rowmax(dl)))
            ang[idx], anl[idx] = dev, dlo
            is_layer[idx] = True
        return phi, ang, phl, anl, is_layer

    def _phi_parts(P, R, ent, trans):
        """横断辺ごとの (φ, φ の下限)。下限 = φ − (辺と入った面の法線の向きの丸めの上限)、0 で打ち切り。"""
        E = P[:, [v for _, v in trans]] - P[:, [u for u, _ in trans]]
        ph = _phi_deg(E, _face_normals(P[:, list(ent)]))
        al = _round_allow_deg(np.linalg.norm(E, axis=2), R[:, [u for u, _ in trans]], R[:, [v for _, v in trans]])
        F = P[:, list(ent)]; RF = R[:, list(ent)]
        fl = np.linalg.norm(np.roll(F, -1, axis=1) - F, axis=2)
        fa = _round_allow_deg(fl, RF, np.roll(RF, -1, axis=1)).max(axis=1)
        return ph, np.where(np.isnan(ph), np.nan, np.maximum(ph - al - fa[:, None], 0.0))

    def _angle_dev_parts(Q, RQ):
        """多角形の角ごとの (内角の 90° からのずれ, その下限)。下限は角を挟む 2 辺の向きの丸めの上限を差し引く。"""
        A = _poly_angles(Q)
        la = np.linalg.norm(np.roll(Q, 1, axis=1) - Q, axis=2)
        lb = np.linalg.norm(np.roll(Q, -1, axis=1) - Q, axis=2)
        al = _round_allow_deg(la, RQ, np.roll(RQ, 1, axis=1)) + _round_allow_deg(lb, RQ, np.roll(RQ, -1, axis=1))
        dv = np.abs(A - 90.0)
        return dv, np.where(np.isnan(A), np.nan, np.maximum(dv - al, 0.0))

    def centroid(c):
        return P3[conn[starts[c]:offs[c]]].mean(axis=0)

    def fmt_xyz(x):
        return "(%.6g, %.6g, %.6g)" % tuple(x)

    def centroids(cells):
        """セルの重心 (節点の平均)。節点数ごとにまとめて計算する。"""
        cells = np.asarray(cells, dtype=np.int64)
        out = np.zeros((len(cells), 3))
        if not len(cells):
            return out
        nn = offs[cells] - starts[cells]
        for m in np.unique(nn):
            i = np.flatnonzero(nn == m)
            out[i] = P3[conn[starts[cells[i]][:, None] + np.arange(int(m))[None, :]]].mean(axis=1)
        return out

    # --- 壁の面 → 隣接セル (physID ごと。対応づけられない physID は判定不能) ---
    resolved, unresolved = {}, {}
    for pid in sorted(walls):
        name, kind, faces = walls[pid]
        if faces is None:
            unresolved[pid] = ("  壁 %s (physID %d, %s): h5 に境界面が無い (BCONDS/%d の vizBfaceNodes・iPlanes) → 判定不能"
                               % (name, pid, kind, pid), "physID %d の境界面が無い" % pid)
            continue
        c0, l0 = lookup(faces)
        bad = int((c0 < 0).sum())
        if bad:
            unresolved[pid] = ("  壁 %s (physID %d, %s): 境界面 %d 個を primal セルに対応づけられない → 判定不能" % (name, pid, kind, bad),
                               "physID %d の面 %d 個が primal セルと対応しない" % (pid, bad))
            continue
        resolved[pid] = (c0, l0)

    # --- 角部 (幾何の角) と境界のラベルの継ぎ目 (#8c) ---
    # 角部 = 節点に接する壁の面の法線どうしの角の最大 (折れ角) が corner-angle 以上の節点。physID は見ない
    # (同一直線上の壁を 2 つの physID に分けた継ぎ目は角部ではない)。法線は隣接セルの重心の側 (領域の内側) へ向けてから
    # 比べる (向きを揃えずに |cos| で比べると、法線がほぼ逆向きに折れる鋭い角を見落とす)。角部は報告を分けるだけで、
    # 緩和の判定・近壁の WARN からは外さない (45° は分類の暫定値で、免除の根拠にしない)。
    inc_node, inc_n = [], []
    wall_count = np.zeros(nnodes, dtype=np.int64)
    for pid, (c0, _) in resolved.items():
        faces = walls[pid][2]
        if not len(faces):
            continue
        wall_count[np.unique(np.concatenate(faces))] += 1
        sz = np.array([len(fc) for fc in faces])
        cc = centroids(c0)
        for k in np.unique(sz):
            i = np.flatnonzero(sz == k)
            F = np.stack([np.asarray(faces[j], dtype=np.int64) for j in i])
            PF = P3[F]
            n = _face_normals(PF)
            ln = np.linalg.norm(n, axis=1)
            sg = np.sign(np.einsum("ij,ij->i", n, cc[i] - PF.mean(axis=1)))
            ok = (ln > 0) & (sg != 0)
            if ok.any():
                nu = n[ok] * (sg[ok] / ln[ok])[:, None]
                inc_node.append(F[ok].ravel()); inc_n.append(np.repeat(nu, int(k), axis=0))
    fold = np.zeros(nnodes)
    if inc_node:
        pn = np.concatenate(inc_node); nv = np.concatenate(inc_n)
        o = np.argsort(pn, kind="stable")
        pn, nv = pn[o], nv[o]
        un, st, cnt = np.unique(pn, return_index=True, return_counts=True)
        for s_ in np.unique(cnt):
            if s_ < 2:
                continue
            g = np.flatnonzero(cnt == s_)
            V = nv[st[g][:, None] + np.arange(int(s_))[None, :]]                 # (節点, 面, 3)
            D = np.einsum("gik,gjk->gij", V, V).reshape(len(g), -1).min(axis=1)
            fold[un[g]] = np.degrees(np.arccos(np.clip(D, -1.0, 1.0)))
    corner_node = fold >= a.corner_angle
    seam_node = (wall_count >= 2) & ~corner_node
    corner_cell = np.maximum.reduceat(corner_node[conn].astype(np.int8), starts).astype(bool) if ncells else np.zeros(0, bool)

    # --- 出力の見出し ---
    summary = {"walls": [], "phi_max": None, "phi_max_regular": None, "phi_max_corner": None, "n_walls": 0}
    say("wall orthogonality (--bcond): 壁 %d 本 (kind %s)、構成の判定 φ <= %.3g°・内角 90±%.3g° (暫定; 精度の限界ではない)、"
        "角部の折れ角 %.3g° (分類の暫定値)" % (len(walls), "/".join(WALL_KINDS), a.phi_construct, a.theta_tol, a.corner_angle))
    say("  φ・内角のずれは格納された座標 (ソルバが使う形) で測る。節点の位置の丸めの上限 %.3g (格納 %s)"
        % ((float(R3.max()) if nnodes else 0.0), np.asarray(coord3).dtype))
    say("  AR の緩和・近壁の WARN は列ごと (壁の面 1 つから壁の外へ向かいの面を渡る層の並び) に判定する。physID は報告の単位")
    if not walls:
        say("  対象の壁 0 (bcond に wall・wall_isothermal が無い)")
    nc_nodes, n_seam = int(corner_node.sum()), int(seam_node.sum())
    kc = int(np.argmax(np.where(corner_node, fold, -1.0))) if nc_nodes else None
    if walls:
        say("  角部 (壁の面の法線の折れ >= %.3g°; 緩和の判定・近壁の WARN からは外さない): 節点 %d%s"
            % (a.corner_angle, nc_nodes, ("、折れ角 最大 %.1f° @ %s" % (float(fold[kc]), fmt_xyz(P3[kc]))) if nc_nodes else ""))
    if n_seam:
        ks = int(np.flatnonzero(seam_node)[0])
        say("  境界のラベルの継ぎ目 (壁の physID が変わるだけで折れ < %.3g°): 節点 %d (例 %s、折れ %.2f°) → 角部として扱わない"
            % (a.corner_angle, n_seam, fmt_xyz(P3[ks]), float(fold[ks])))
    summary["corners"] = {"n_nodes": nc_nodes, "fold_max": None if kc is None else float(fold[kc]),
                          "at": None if kc is None else [float(x) for x in P3[kc]], "n_seam_nodes": n_seam,
                          "corner_angle": a.corner_angle}

    # --- 壁ごと (報告の単位): 層 1 の横断辺の φ と位置 ---
    all_phi = {"通常": None, "角部": None}
    for pid in sorted(walls):
        name, kind, faces = walls[pid]
        if pid in unresolved:
            say(unresolved[pid][0])
            indet.append(unresolved[pid][1])
            continue
        c0, l0 = resolved[pid]
        _, _, _, _, lay0 = layer_metrics(c0, l0)
        edge_phi, edge_lo, edge_xyz, edge_corner = [], [], [], []
        code = types[c0] * 16 + l0
        for cval in np.unique(code):
            t, lf = int(cval // 16), int(cval % 16)
            g = _layer_geometry(t, lf, topo3d)
            if g is None:
                continue
            ent, trans, _ = g
            idx = np.flatnonzero(code == cval)
            nn = int(offs[c0[idx[0]]] - starts[c0[idx[0]]])
            ids = conn[starts[c0[idx]][:, None] + np.arange(nn)[None, :]]
            P = P3[ids]
            ph, pl = _phi_parts(P, R3[ids], ent, trans)
            for j, (u, _) in enumerate(trans):
                edge_phi.append(ph[:, j]); edge_lo.append(pl[:, j]); edge_xyz.append(P[:, u])
                edge_corner.append(corner_cell[c0[idx]])
        n_out = int((~lay0).sum())
        out_types = {}
        for t in types[c0[~lay0]]:
            out_types[VTK_NAME.get(int(t), str(int(t)))] = out_types.get(VTK_NAME.get(int(t), str(int(t))), 0) + 1
        rec = {"physID": pid, "name": name, "kind": kind, "faces": len(faces), "target_faces": int(lay0.sum()),
               "excluded_faces": n_out, "excluded_types": out_types}
        say("  壁 %s (physID %d, %s): 境界面 %d、対象 %d、対象外 %d%s" % (
            name, pid, kind, len(faces), int(lay0.sum()), n_out,
            (" (" + ", ".join("%s %d" % kv for kv in sorted(out_types.items())) + ")") if n_out else ""))
        if edge_phi:
            ephi = np.concatenate(edge_phi); elo = np.concatenate(edge_lo)
            exyz = np.concatenate(edge_xyz); ecor = np.concatenate(edge_corner)
            for lab, m in (("通常", ~ecor), ("角部", ecor)):
                if not m.any():
                    continue
                v = ephi[m]
                if np.isnan(v).any():
                    say("    %s: 対象辺 %d のうち退化 (長さ 0) %d → 判定不能" % (lab, int(m.sum()), int(np.isnan(v).sum())))
                    indet.append("physID %d の退化した辺" % pid)
                    v = np.nan_to_num(v, nan=0.0)
                k = int(np.argmax(v))
                rmax = float(np.nanmax(v - np.nan_to_num(elo[m], nan=0.0)))      # 丸めの上限 (φ − 下限) の最大
                bins = "  ".join(">%g°: %d" % (b, int((v > b).sum())) for b in PHI_BINS)
                say("    %s: 対象辺 %d (面ごと)、φ 最大 %.2f° @ %s  %s  (丸めの上限 最大 %.2f°)"
                    % (lab, int(m.sum()), float(v[k]), fmt_xyz(exyz[m][k]), bins, rmax))
                rec["edges_" + ("corner" if lab == "角部" else "regular")] = {
                    "n": int(m.sum()), "phi_max": float(v[k]), "at": [float(x) for x in exyz[m][k]],
                    "bins": {str(b): int((v > b).sum()) for b in PHI_BINS}, "rounding_max_deg": rmax}
                all_phi[lab] = float(v[k]) if all_phi[lab] is None else max(all_phi[lab], float(v[k]))
        summary["walls"].append(rec)
    summary["n_walls"] = len(walls)
    summary["phi_max_regular"], summary["phi_max_corner"] = all_phi["通常"], all_phi["角部"]
    have = [x for x in all_phi.values() if x is not None]
    summary["phi_max"] = max(have) if have else None              # 角部を含む (角部を最大値から外さない)

    # --- 列ごとの行進 (#8c): 壁の面 1 つ = 列 1 本。入った面の向かいの面を渡って壁の外へ進む (層 1, 2, ...) ---
    # セルは最も近い (層の番号が最小の) 列に属する。同じ層の番号で 2 本以上の列が届くときは全部の列で条件を満たすこと。
    # 列は層でない型 (三角形など) か境界で終わる。境界から始まる列は同じ (セル, 入った面) に戻らないので必ず終わる。
    tmask = (ars > AR_BASE) & (ars <= a.ar_max) if a.ar_max > AR_BASE else np.zeros(ncells, dtype=bool)
    minD = np.zeros(ncells, dtype=np.int64)                      # 0 = どの列も届かない
    VER = np.zeros(ncells, dtype=bool); VLO = np.zeros(ncells, dtype=bool); BADc = np.zeros(ncells, dtype=bool)
    PH1 = np.full(ncells, np.nan); PL1 = np.full(ncells, np.nan); LY1 = np.zeros(ncells, dtype=bool)
    tp = {k: [] for k in ("c", "pid", "phi", "phi_at", "ang", "ang_at", "depth")}
    if resolved:
        pids_ok = sorted(resolved)
        C0 = np.concatenate([resolved[p][0] for p in pids_ok])
        L0 = np.concatenate([resolved[p][1] for p in pids_ok])
        PIDc = np.concatenate([np.full(len(resolved[p][0]), p, dtype=np.int64) for p in pids_ok])
        _, first = np.unique(C0 * 16 + L0, return_index=True)    # 同じ面が 2 つの physID に入っていても列は 1 本
        first = np.sort(first)
        C0, L0, PIDc = C0[first], L0[first], PIDc[first]
        cur_col = np.arange(len(C0)); cur_c = C0.copy(); cur_l = L0.copy()
        pok = np.ones(len(C0), dtype=bool); plo = np.ones(len(C0), dtype=bool); pbad = np.zeros(len(C0), dtype=bool)
        pphi = np.full(len(C0), -np.inf); pphi_at = np.full(len(C0), -1, dtype=np.int64)
        pang = np.full(len(C0), -np.inf); pang_at = np.full(len(C0), -1, dtype=np.int64)
        n_left = int(tmask.sum())
        depth = 1
        while len(cur_c) and depth <= ncells * maxlf:
            ph, an, pl, al, ly = layer_metrics(cur_c, cur_l)
            fin = np.isfinite(ph) & np.isfinite(an)
            with np.errstate(invalid="ignore"):
                pok &= ly & fin & (ph <= a.phi_construct) & (an <= a.theta_tol)
                plo &= ly & np.isfinite(pl) & np.isfinite(al) & (pl <= a.phi_construct) & (al <= a.theta_tol)
            pbad |= ~(ly & fin)
            phv = np.where(np.isfinite(ph), ph, -np.inf); anv = np.where(np.isfinite(an), an, -np.inf)
            u = phv > pphi
            pphi = np.where(u, phv, pphi); pphi_at = np.where(u, cur_c, pphi_at)
            u = anv > pang
            pang = np.where(u, anv, pang); pang_at = np.where(u, cur_c, pang_at)
            new = minD[cur_c] == 0
            if new.any():
                cn = cur_c[new]
                minD[cn] = depth
                VER[cn] = True; np.logical_and.at(VER, cn, pok[new])
                VLO[cn] = True; np.logical_and.at(VLO, cn, plo[new])
                np.logical_or.at(BADc, cn, pbad[new])
                if depth == 1:
                    np.fmax.at(PH1, cn, ph[new]); np.fmax.at(PL1, cn, pl[new]); np.logical_or.at(LY1, cn, ly[new])
                tm = tmask[cn]
                if tm.any():
                    sel = np.flatnonzero(new)[tm]
                    tp["c"].append(cur_c[sel]); tp["pid"].append(PIDc[cur_col[sel]])
                    tp["phi"].append(pphi[sel]); tp["phi_at"].append(pphi_at[sel])
                    tp["ang"].append(pang[sel]); tp["ang_at"].append(pang_at[sel])
                    tp["depth"].append(np.full(len(sel), depth, dtype=np.int64))
                    n_left -= len(np.unique(cn[tm]))
            if n_left <= 0:
                break                                            # 層 1 と全部の対象のセルが済んだ (この先は結果を変えない)
            op = opp_tab[types[cur_c], cur_l]
            ok = ly & (op >= 0)
            nc = np.where(ok, nbr_c[cur_c, np.maximum(op, 0)], -1)
            nl = np.where(ok, nbr_l[cur_c, np.maximum(op, 0)], -1)
            m = nc >= 0
            cur_col, cur_c, cur_l = cur_col[m], nc[m], nl[m]
            pok, plo, pbad, pphi, pphi_at, pang, pang_at = pok[m], plo[m], pbad[m], pphi[m], pphi_at[m], pang[m], pang_at[m]
            depth += 1
    tp = {k: (np.concatenate(v) if v else np.zeros(0, dtype=(np.float64 if k in ("phi", "ang") else np.int64)))
          for k, v in tp.items()}

    # --- AR の緩和 (1000 < AR <= --ar-max): 壁の面から対象のセルまでの同じ列の全層 ---
    n_unverified = 0
    relax = {"n_target": 0, "n_verified": 0, "n_unverified": 0}
    verified_mask = np.zeros(ncells, dtype=bool)
    if a.ar_max > AR_BASE:
        T = np.flatnonzero(tmask)
        reached = minD[T] > 0
        ver = reached & VER[T]
        verified_mask[T[ver]] = True
        why = np.full(len(T), "", dtype=object)
        why[~reached] = "層の外"
        r = reached & ~ver
        why[r & BADc[T]] = "層でない型・退化"
        why[r & ~BADc[T] & VLO[T]] = "層の条件 (丸めの範囲)"
        why[r & ~BADc[T] & ~VLO[T]] = "層の条件"
        reasons = {}
        for k in ("層の外", "層でない型・退化", "層の条件", "層の条件 (丸めの範囲)"):
            if int((why == k).sum()):
                reasons[k] = int((why == k).sum())
        n_unverified = int((~ver).sum())
        relax.update(n_target=int(len(T)), n_verified=int(ver.sum()), n_unverified=n_unverified, reasons=reasons)
        per_pid = []
        for pid in np.unique(tp["pid"]):
            m = tp["pid"] == pid
            cs = np.unique(tp["c"][m])
            nv_ = int(VER[cs].sum())
            nr_ = int((~VER[cs] & ~BADc[cs] & VLO[cs]).sum())
            d = tp["depth"][m]
            kp = int(np.argmax(tp["phi"][m])); ka = int(np.argmax(tp["ang"][m]))
            fp, fa = float(tp["phi"][m][kp]), float(tp["ang"][m][ka])
            name = walls[int(pid)][0]
            say("  AR 緩和の列 (壁 %s physID %d): 対象 %d セル (層 %d..%d)、壁から対象までの全層の φ 最大 %s、内角のずれ最大 %s"
                " → 確認 %d、確認できない %d%s"
                % (name, int(pid), len(cs), int(d.min()), int(d.max()),
                   ("%.2f° @ %s" % (fp, fmt_xyz(centroid(int(tp["phi_at"][m][kp]))))) if np.isfinite(fp) else "— (層でない型・退化)",
                   ("%.2f° @ %s" % (fa, fmt_xyz(centroid(int(tp["ang_at"][m][ka]))))) if np.isfinite(fa) else "—",
                   nv_, len(cs) - nv_,
                   (" (うち %d セルは丸めの上限を差し引けば満たす: 格納精度で確かめられない。FP64 の変換か層を厚く)" % nr_) if nr_ else ""))
            per_pid.append({"physID": int(pid), "n_target": int(len(cs)), "n_verified": nv_, "depth": [int(d.min()), int(d.max())],
                            "phi_max": fp if np.isfinite(fp) else None, "ang_dev_max": fa if np.isfinite(fa) else None,
                            "n_within_rounding": nr_})
        relax["per_wall"] = per_pid
        say("AR 緩和 (%.0f < AR <= %.0f、列ごとに壁から対象のセルまでの全層): 対象 %d セル → 壁法線の構造層として確認 %d、確認できない %d%s"
            % (AR_BASE, a.ar_max, len(T), int(ver.sum()), n_unverified,
               (" (" + ", ".join("%s %d" % kv for kv in reasons.items()) + ")") if n_unverified else ""))
        if n_unverified:
            bad = T[~ver]
            k = bad[int(np.argmax(ars[bad]))]
            say("  確認できないセルの AR 最大 %.1f @ %s (件数によらず FAIL: SOFT-PASS を当てない)" % (float(ars[k]), fmt_xyz(centroid(k))))
    summary["relax"] = relax

    # --- 近壁 (層 1) の φ の WARN (AR <= 1000)。角部のセルも含める (角部は報告を分けるだけ) ---
    with np.errstate(invalid="ignore"):
        wmask = (minD == 1) & LY1 & (ars <= AR_BASE) & (np.nan_to_num(PH1, nan=0.0) > a.phi_construct)
    n_warn = int(wmask.sum())
    n_warn_corner = int((wmask & corner_cell).sum())
    n_warn_round = int((wmask & (np.nan_to_num(PL1, nan=0.0) <= a.phi_construct)).sum())
    wmax = wat = cmax = cat = None
    if n_warn:
        cs = np.flatnonzero(wmask)
        k = cs[int(np.argmax(PH1[cs]))]
        wmax, wat = float(PH1[k]), centroid(k)
        if n_warn_corner:
            cs = np.flatnonzero(wmask & corner_cell)
            k = cs[int(np.argmax(PH1[cs]))]
            cmax, cat = float(PH1[k]), centroid(k)
        say("WARN 近壁 (層 1、AR <= %.0f) で φ > %.3g°: %d セル、φ 最大 %.2f° @ %s%s%s"
            % (AR_BASE, a.phi_construct, n_warn, wmax, fmt_xyz(wat),
               ("  [うち角部 %d セル、φ 最大 %.2f° @ %s]" % (n_warn_corner, cmax, fmt_xyz(cat))) if n_warn_corner else "",
               ("  [うち %d セルは丸めの上限の範囲]" % n_warn_round) if n_warn_round else ""))
    summary["warn_near_wall"] = {"n": n_warn, "n_corner": n_warn_corner, "n_within_rounding": n_warn_round, "phi_max": wmax,
                                 "at": None if wat is None else [float(x) for x in wat],
                                 "corner_phi_max": cmax, "corner_at": None if cat is None else [float(x) for x in cat]}
    masks = {"target": tmask, "verified": verified_mask, "warn": wmask, "corner_cell": corner_cell, "layer": minD}
    return {"lines": lines, "n_unverified": n_unverified, "warn": n_warn, "indeterminate": indet, "summary": summary,
            "masks": masks}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--ar-max", type=float, default=1000.0,
                    help="AR の上限 (既定 1000、最大 %.0f)。1000 を超える値は --bcond があれば壁法線の構造層のセルだけに効く" % AR_RELAX_CAP)
    ap.add_argument("--skew-max", type=float, default=0.9)
    ap.add_argument("--mode", choices=["auto", "2d", "3d"], default="auto")
    ap.add_argument("--loop", action="store_true", help="旧のセル毎ループで計算 (検証用, 遅い)")
    ap.add_argument("--no-vizmesh", action="store_true", help="/VIZMESH があっても /MESH/CONNE を測る")
    ap.add_argument("--bcond", default=None,
                    help="bcondConfig.yaml。壁 (wall・wall_isothermal) の physID を h5 の BCONDS/<physID> と対応づけ、"
                         "壁直交性 φ と AR の緩和の条件を判定する")
    ap.add_argument("--phi-construct", type=float, default=PHI_CONSTRUCT_DEFAULT,
                    help="構成の判定の φ の上限 [deg] (暫定 %.3g)。層が壁法線に沿って作られたかを分ける値で、"
                         "精度の限界ではない (どこまで傾いてよいかは壁応力・熱流束の A/B で別に決める)" % PHI_CONSTRUCT_DEFAULT)
    ap.add_argument("--theta-tol", type=float, default=THETA_TOL_DEFAULT,
                    help="緩和の対象の層の内角の 90° からのずれの上限 [deg] (暫定 %.3g、構成の判定)" % THETA_TOL_DEFAULT)
    ap.add_argument("--corner-angle", type=float, default=CORNER_ANGLE_DEFAULT,
                    help="壁の面の法線がこれ以上折れる節点を角部として分けて数える [deg] (既定 %.3g)" % CORNER_ANGLE_DEFAULT)
    ap.add_argument("--no-cache", action="store_true", help="<h5>.quality.json の結果を再利用しない (書きもしない)")
    ap.add_argument("--json", default=None, help="判定の要約を JSON で書く (--bcond のとき)")
    a = ap.parse_args()
    if a.ar_max > AR_RELAX_CAP:
        ap.error("--ar-max は %.0f まで (AGENTS.md の AR 緩和の上限)。指定 %g" % (AR_RELAX_CAP, a.ar_max))

    out_lines = []

    def say(msg):
        print(msg)
        out_lines.append(msg)

    walls, bc_kinds, cache_key, cache_path = None, None, None, None
    with h5py.File(a.h5, "r") as f:
        coord3 = np.array(f["MESH/COORD"]).reshape(-1, 3)
        if "VIZMESH/CONNE" in f and not a.no_vizmesh:
            # node (median-dual) 変換: primal セル (可視化トポロジ) を測る
            conne = np.array(f["VIZMESH/CONNE"]); ncells = int(f["VIZMESH"].attrs["nVizCells"])
            src_line = f"source: /VIZMESH/CONNE (primal cells of a node-converted mesh, {ncells} cells)"
        else:
            conne = np.array(f["MESH/CONNE"])
            ncells = f["VALUE/ro"].shape[0] if "VALUE/ro" in f else None
            src_line = None
        if a.bcond is not None:
            try:
                bc_kinds = read_bcond_kinds(a.bcond)
            except Exception as e:                       # noqa: BLE001
                print(f"--bcond {a.bcond} を読めない: {e}")
                print("VERDICT: FAIL (判定不能: 壁の bcond が分からない)")
                sys.exit(1)
            wpids = sorted(p for p, (_, k) in bc_kinds.items() if k in WALL_KINDS)
            faces = read_boundary_faces(f, wpids)
            walls = {p: (bc_kinds[p][0], bc_kinds[p][1], faces[p]) for p in wpids}
            if not a.no_cache:
                cache_path = a.h5 + ".quality.json"
                cache_key = {
                    "checker": CHECKER_VERSION, "checker_sha256": checker_sha(),
                    "mesh_invariant_sha256": mesh_invariant_sha(f),
                    "walls": {str(p): (None if faces[p] is None else _faces_sha(faces[p])) for p in wpids},
                    "args": {"ar_max": a.ar_max, "skew_max": a.skew_max, "mode": a.mode, "no_vizmesh": a.no_vizmesh,
                             "phi_construct": a.phi_construct, "theta_tol": a.theta_tol, "corner_angle": a.corner_angle},
                }
            h5_kinds = {}
            if "BCONDS" in f:
                for k in f["BCONDS"]:
                    kk = f["BCONDS"][k].attrs.get("bcondKind")
                    h5_kinds[int(k)] = kk.decode() if isinstance(kk, bytes) else (None if kk is None else str(kk))

    def finish(code, verdict, summary=None):
        """--bcond のとき: 結果を <h5>.quality.json に書き (鍵つき)、--json を書いて終わる。"""
        if cache_path is not None:
            tmp = f"{cache_path}.tmp{os.getpid()}"
            try:
                with open(tmp, "w") as fh:          # 書きかけを別の検査が読まないよう、別名に書いてから置き換える
                    json.dump({"key": cache_key, "result": {"exit": code, "verdict": verdict, "lines": out_lines,
                                                            "summary": summary}}, fh, ensure_ascii=False, indent=1)
                os.replace(tmp, cache_path)
            except OSError as e:
                print(f"(結果を {cache_path} に書けない: {e})")
        if a.json:
            with open(a.json, "w") as fh:
                json.dump({"verdict": verdict, "exit": code, "reused": False, "summary": summary}, fh, ensure_ascii=False, indent=1)
        sys.exit(code)

    if cache_key is not None:
        prev, why = None, None
        if os.path.exists(cache_path):
            try:
                with open(cache_path) as fh:
                    prev = json.load(fh)
            except Exception as e:                       # noqa: BLE001
                why = f"読めない ({e})"
        else:
            why = "記録なし"
        if prev is not None:
            pk = prev.get("key") or {}
            diff = [k for k in ("checker", "checker_sha256", "mesh_invariant_sha256", "walls", "args") if pk.get(k) != cache_key[k]]
            if not diff and "result" in prev:
                res = prev["result"]
                print(f"結果の再利用: {cache_path} (鍵が一致: メッシュの不変部分・壁の physID と面・閾値・検査器の版)")
                for ln in res.get("lines", []):
                    print(ln)
                if a.json:
                    with open(a.json, "w") as fh:
                        json.dump({"verdict": res.get("verdict"), "exit": res.get("exit"), "reused": True,
                                   "summary": res.get("summary")}, fh, ensure_ascii=False, indent=1)
                sys.exit(int(res.get("exit", 1)))
            why = "鍵が違う: " + ", ".join(diff)
        print(f"結果の再利用: なし ({why})")

    if src_line is not None:
        say(src_line)
    if ncells is None:
        # count from CONNE
        from res_h5_to_vtu import _count_cells
        ncells = _count_cells(conne)
    conn, offs, types = parse_conne(conne, ncells)

    # モード判定: 3D セル型を含み z 平面が 3 枚以上なら 3D (x-y 射影が破綻するため)
    if a.mode == "auto":
        nz = len(np.unique(np.round(coord3[:, 2], 6)))
        is3d = bool(np.isin(types, list(VTK_3D)).any()) and nz > 2
    else:
        is3d = (a.mode == "3d")
    say(f"mode: {'3D (per-face skew / edge AR)' if is3d else '2D (x-y projection)'}")

    coord = coord3 if is3d else coord3[:, :2]
    if a.loop:
        ars = np.empty(ncells); skews = np.empty(ncells)
        s = 0
        for i, o in enumerate(offs):
            pts = coord[conn[s:o]]; s = o
            ars[i], skews[i] = (cell_metrics_3d(pts, types[i]) if is3d else cell_metrics(pts))
    else:
        ars, skews = metrics_vectorized(coord, np.asarray(conn), np.asarray(offs), np.asarray(types), is3d)

    # --- メッシュの**成立性**を、形状品質の外れ値許容より前に無条件で検査する ---
    # (2026-09-19 codex: NaN 座標が 1 セルあっても `SOFT-PASS` exit 0、体積ゼロの四面体は
    #  AR 1.732 / skew 0.500 で通常合格していた。形状の良し悪し以前の問題は許容しない。)
    fatal = []
    nonfinite = int(np.count_nonzero(~np.isfinite(coord3)))
    if nonfinite:
        bad_nodes = int(np.count_nonzero(~np.isfinite(coord3).all(axis=1)))
        fatal.append("非有限座標: 成分 %d 個 / 節点 %d 個" % (nonfinite, bad_nodes))
    if is3d:
        vols, n_inverted = cell_volumes(coord3, np.asarray(conn), np.asarray(offs), np.asarray(types))
        if n_inverted:
            fatal.append("六面体の頂点 Jacobian の符号が混在 (局所反転) %d 個" % n_inverted)
        nv_bad = int(np.count_nonzero(~np.isfinite(vols)))
        nz_bad = int(np.count_nonzero(np.abs(vols) <= 0.0))
        neg = int(np.count_nonzero(vols < 0.0))
        pos = int(np.count_nonzero(vols > 0.0))
        if nv_bad:
            fatal.append("体積が非有限のセル %d 個" % nv_bad)
        if nz_bad:
            fatal.append("体積ゼロのセル %d 個" % nz_bad)
        if neg and pos:
            fatal.append("体積の符号が混在 (正 %d / 負 %d) = 向きの不整合" % (pos, neg))
        good = vols[np.isfinite(vols) & (vols != 0.0)]
        if good.size:
            say("cell volume  : |V| min=%.3e  median=%.3e  max=%.3e  (0 or 非有限: %d)"
                % (np.abs(good).min(), np.median(np.abs(good)), np.abs(good).max(),
                   nv_bad + nz_bad))
    if fatal:
        say("FATAL: " + " / ".join(fatal))
        say("VERDICT: FAIL (メッシュとして成立していない。形状品質の外れ値許容は適用しない)")
        if a.bcond is not None:
            finish(1, "FAIL", {"fatal": fatal})
        sys.exit(1)

    ar_bad = int((ars > a.ar_max).sum())
    sk_bad = int((skews > a.skew_max).sum())
    say(f"cells: {ncells}")
    say(f"aspect ratio : max={ars.max():.1f}  p99={np.percentile(ars,99):.1f}  "
        f"mean={ars.mean():.1f}  >{a.ar_max:.0f}: {ar_bad} ({100*ar_bad/ncells:.2f}%)")
    say(f"skewness     : max={skews.max():.3f}  p99={np.percentile(skews,99):.3f}  "
        f"mean={skews.mean():.3f}  >{a.skew_max:.2f}: {sk_bad} ({100*sk_bad/ncells:.2f}%)")
    if a.bcond is None:
        if a.ar_max > AR_BASE:
            say(f"注意: --bcond なしの --ar-max {a.ar_max:.0f} は旧判定 (緩和の条件 = 壁法線の構造層を確かめずに "
                f"AR <= {a.ar_max:.0f} を通す)。関門・新しいメッシュは --bcond bcondConfig.yaml を渡すこと")
        ok = (ars.max() <= a.ar_max) and (skews.max() <= a.skew_max)
        # 少数の外れ値は許容するソフト判定も併記
        soft = (ar_bad <= 0.001 * ncells) and (sk_bad <= 0.001 * ncells)
        say("注意: --bcond なしの判定は投入の関門の合格の根拠にならない (壁直交性・AR の緩和の条件を確かめていない。"
            "関門は --bcond bcondConfig.yaml で判定する)")
        say(f"VERDICT: {'PASS' if ok else ('SOFT-PASS (<0.1% outliers)' if soft else 'FAIL')} "
            f"(AR<= {a.ar_max:.0f}, skew<= {a.skew_max:.2f})")
        sys.exit(0 if (ok or soft) else 1)

    # --- 壁直交性と AR の緩和の条件 (--bcond) ---
    mism = sorted(p for p in walls if p in h5_kinds and h5_kinds[p] not in WALL_KINDS)
    mism += sorted(p for p, (_, k) in bc_kinds.items() if k not in WALL_KINDS and h5_kinds.get(p) in WALL_KINDS)
    if mism:
        say("注意: 変換時の bcondKind と --bcond の kind の壁・非壁が違う physID: %s (メッシュの来歴の照合は関門の別の検査)"
            % ", ".join("%d (h5 %s / bcond %s)" % (p, h5_kinds.get(p), bc_kinds[p][1]) for p in sorted(set(mism))))
    missing = sorted(set(bc_kinds) - set(h5_kinds)) if h5_kinds else []
    if missing:
        say("注意: --bcond の physID が h5 の BCONDS に無い: %s" % ", ".join(str(p) for p in missing))
    wo = wall_orthogonality(coord3, conn, offs, types, walls, ars, a, rnode=node_rounding_radius(coord3))
    for ln in wo["lines"]:
        say(ln)
    n_unv = wo["n_unverified"]
    # AR > --ar-max (AR > 5000 を含む) は件数によらず FAIL (#8c: --bcond では AR の外れ値に SOFT-PASS を当てない。
    # 0.1 % 未満の 1 セルの AR 6000.3 が SOFT-PASS・終了コード 0 で通った: codex diagnose 2026-10-14)
    n_cap = int((ars > AR_RELAX_CAP).sum())
    if ar_bad:
        k = int(np.argmax(ars))
        s0 = int(offs[k - 1]) if k else 0
        ctr = coord3[np.asarray(conn)[s0:int(offs[k])]].mean(axis=0)
        say("AR > %.0f (--ar-max) のセル %d%s、AR 最大 %.1f @ (%.6g, %.6g, %.6g) → 件数によらず FAIL "
            "(--bcond では AR の外れ値に SOFT-PASS を当てない%s)"
            % (a.ar_max, ar_bad, (" (うち AR > %.0f: %d)" % (AR_RELAX_CAP, n_cap)) if a.ar_max < AR_RELAX_CAP else "",
               float(ars[k]), *ctr,
               "" if a.ar_max > AR_BASE else "。壁法線の構造層を AR 5000 まで認めるには --ar-max を上げて緩和の条件を確かめる"))
    ok = (ar_bad == 0) and (sk_bad == 0) and (n_unv == 0)
    soft = (n_unv == 0) and (ar_bad == 0) and (sk_bad <= 0.001 * ncells)      # SOFT-PASS は skew の外れ値だけ
    verdict = "PASS" if ok else ("SOFT-PASS" if soft else "FAIL")
    if verdict != "FAIL" and wo["indeterminate"]:
        verdict = "FAIL"
        say("判定不能: " + "; ".join(wo["indeterminate"]))
    if verdict != "FAIL" and wo["warn"]:
        verdict = "WARN"
    detail = {"PASS": "", "SOFT-PASS": " (<0.1% skew outliers)", "FAIL": "",
              "WARN": " (近壁の φ > %.3g°%s)" % (a.phi_construct, "; <0.1% skew outliers" if not ok else "")}[verdict]
    relax_txt = (f", {AR_BASE:.0f} < AR <= {a.ar_max:.0f} は壁法線の構造層だけ" if a.ar_max > AR_BASE else "")
    say(f"VERDICT: {verdict}{detail} (AR<= {min(a.ar_max, AR_BASE):.0f}{relax_txt}, skew<= {a.skew_max:.2f}, "
        f"構成 φ<= {a.phi_construct:.3g}°)")
    summary = dict(wo["summary"], ar_bad=ar_bad, ar_over_cap=n_cap, skew_bad=sk_bad, ncells=int(ncells),
                   indeterminate=wo["indeterminate"])
    finish(0 if verdict != "FAIL" else 1, verdict, summary)


if __name__ == "__main__":
    main()
