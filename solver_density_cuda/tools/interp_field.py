#!/usr/bin/env python3
"""過去 h5 (任意メッシュ) の場を最近傍interpolateで新メッシュ入力 h5 に貼る (cross-mesh restart)。

メッシュを変えた (quad↔tri↔構造, 解像度変更) ときに、収束済みの場を初期値として移植する
標準ツール。**メッシュ変更後は uniform 初期値から始めず本ツールで restart する**
(超音速/衝撃波/SST は uniform IC から発散しやすい)。同一メッシュの restart は restart_field.py。

- SRC: res_*.h5 (primitives P,T,Ux,Uy,Uz,ro[,k,omega]) でも input h5 (conserved) でも可。
- DST: convertGmshToForge 直後の新メッシュ input h5。
- 各 DST セル重心に最も近い SRC セル重心の値を採用 (scipy cKDTree, 最近傍)。
- 移植する /VALUE/: 保存量 ro,roUx,roUy,roUz,roe・乱流 roK,roOmega・スカラー輸送 roY*・凝縮モーメント rog_*/roQ*_* (あれば)。
- **wall_dist は移植しない** (新メッシュで convert 時に計算済みの値を使う)。

- **化学種の照合**: SRC/DST の隣に `solverConfig.yaml` があれば `physProp.species` (名前と順序) を比較し、
  違えば拒否する (別順序・別種集合の場を黙って index で貼ると組成が入れ替わる)。種を変える restart は
  `tools/convert_species_field.py` (擬似種の展開・名前で移す) を使う。`--force-species` で照合を無視できる。

usage: interp_field.py SRC.h5 DST_input.h5 [--gamma 1.4] [--force-species]
"""
import argparse, sys, os
import numpy as np, h5py
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from res_h5_to_vtu import parse_conne


def centroids(f):
    # DOF の代表座標。node (median-dual) では値数=節点数なので MESH/COORD (節点座標)
    # を最優先する — CELLS/centCoords (双対 CV 重心) は高 AR の μm 級壁 CV で節点から
    # ±100 μm 級に外れ (壁 CV の重心が域外に出る例あり)、最近傍照会に使うと近壁 IC が
    # src の壁値/内部値を交互に拾う市松になる (case/40 y+1 node で発散種になった実害)。
    coord = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, :2]
    if "VALUE/ro" in f and f["VALUE/ro"].shape[0] == coord.shape[0]:
        return coord
    # cell: 入力 h5 は /CELLS/centCoords を持つのでそれを優先 (res_*.h5 は CONNE 経路)
    if "CELLS/centCoords" in f:
        return np.array(f["CELLS/centCoords"]).reshape(-1, 3)[:, :2]
    nc = f["VALUE/ro"].shape[0]
    if nc == coord.shape[0]:
        # node-centered res (median-dual): 値の位置はノード座標そのもの
        # (MESH/CONNE は可視化用 primal トポロジで parse できない)
        return coord
    conn, offs, _ = parse_conne(np.array(f["MESH/CONNE"]), nc)
    c = np.zeros((nc, 2)); s = 0
    for i, o in enumerate(offs):
        c[i] = coord[conn[s:o]].mean(axis=0); s = o
    return c


def _species_names_near(h5path):
    """h5 と同じディレクトリの solverConfig.yaml から physProp.species を返す (無ければ None)。"""
    try:
        from forge_species import species_info
        d = os.path.dirname(os.path.abspath(h5path))
        if not os.path.exists(os.path.join(d, "solverConfig.yaml")):
            return None
        info = species_info(d)
        return list(info["names"]) if info["thermalMethod"] == 2 else None
    except Exception as e:   # noqa: BLE001 — 照合は補助なので読めなければ無視
        print(f"[interp_field] species check skipped for {h5path}: {e}")
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--force-species", action="store_true",
                    help="SRC/DST の physProp.species が違っても index で貼る (通常は convert_species_field.py を使う)")
    a = ap.parse_args(); g = a.gamma

    # 化学種の名前照合 (両側に solverConfig.yaml があるときだけ)。
    src_sp = _species_names_near(a.src); dst_sp = _species_names_near(a.dst)
    if src_sp is not None and dst_sp is not None and src_sp != dst_sp:
        msg = (f"physProp.species が違う: SRC {src_sp} vs DST {dst_sp}. 種の順序/集合が違う場は index コピーできない。"
               f" tools/convert_species_field.py SRC_res.h5 DST_input.h5 --meta DST/species_meta.yaml で名前により移す"
               f" (どうしても index で貼るなら --force-species)。")
        if not a.force_species:
            raise SystemExit("[interp_field] REFUSED: " + msg)
        print("[interp_field] WARNING (--force-species): " + msg)

    with h5py.File(a.src, "r") as s:
        cs = centroids(s); V = s["VALUE"]
        if "P" in V and "Ux" in V:            # res (primitives)
            ro = np.array(V["ro"]); P = np.array(V["P"])
            Ux = np.array(V["Ux"]); Uy = np.array(V["Uy"]); Uz = np.array(V["Uz"])
            # roe: res に保存量 /VALUE/roe があればそれを**そのまま**使う (2026-08-17)。
            # CPG 再構成 P/(γ-1)+½ρu² は thermally-perfect (thermalMethod 2, NASA-9 の絶対
            # 基準 e(T)) では別物になり、移植後に T が数千 K へ跳ぶ (case/42 で実測)。
            # 保存量が無い旧 res のみ CPG 式へフォールバック。
            if "roe" in V:
                roe = np.array(V["roe"])
            else:
                roe = P/(g-1.0) + 0.5*ro*(Ux**2+Uy**2+Uz**2)
            fields = {"ro": ro, "roUx": ro*Ux, "roUy": ro*Uy, "roUz": ro*Uz, "roe": roe}
            if "k" in V and "omega" in V:
                fields["roK"] = ro*np.array(V["k"]); fields["roOmega"] = ro*np.array(V["omega"])
            for key in V:                      # scalar transport Y* -> roY*
                if key.startswith("Y") and key[1:].isdigit():
                    fields["ro"+key] = ro*np.array(V[key])
                if key == "roXi":              # 受動トレーサ (physProp.tracer)
                    fields["roXi"] = np.array(V[key])
                # 凝縮モーメント (原始 g_<s>,Q0_<s>.. → 保存 rog_<s>,roQ0_<s>..)。2026-08-17
                if key.startswith(("g_", "Q0_", "Q1_", "Q2_")):
                    fields["ro"+key] = ro*np.array(V[key])
        else:                                  # input (conserved)
            fields = {n: np.array(V[n]) for n in
                      ["ro","roUx","roUy","roUz","roe","roK","roOmega"] if n in V}
            for key in V:
                if key.startswith("roY") and key[3:].isdigit():
                    fields[key] = np.array(V[key])
                if key.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or key == "roXi":
                    fields[key] = np.array(V[key])

    tree = cKDTree(cs)
    with h5py.File(a.dst, "r+") as d:
        cd = centroids(d)
        _, idx = tree.query(cd)
        moved = []
        for name, arr in fields.items():
            ds = "VALUE/"+name
            if ds in d and name != "wall_dist":
                d[ds][...] = arr[idx].astype(d[ds].dtype); moved.append(name)
            elif name.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or (name.startswith("roY") and name[3:].isdigit()) or name == "roXi":
                # 凝縮モーメントと化学種は convert 直後の入力 h5 に無いので新規作成する (無ければ forge は第 1 種以外を 0 に
                # 初期化し、carrier では rog<=roY_w のクランプで液相が消える: codex 2026-09-16 result M2)。
                # forge は VALUE/<consName> が存在すれば読む (無ければ 0 = dry restart)。2026-08-18 / 2026-09-16
                d.create_dataset(ds, data=arr[idx].astype(d["VALUE/ro"].dtype)); moved.append(name+"(new)")
        print(f"interp {a.src} -> {a.dst}: {len(cd)} dst cells, moved {moved} (wall_dist kept)")


if __name__ == "__main__":
    main()
