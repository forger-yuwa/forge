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

- **化学種の照合 (既定で必須)**: SRC/DST それぞれの隣の `solverConfig.yaml` + `species_db.yaml` (+ `species_meta.yaml`) から
  署名 (`forge_species.species_signature`: 種名と順序・種ごとの MW (rtol 1e-9)・`thermoHrefTemp`・meta の sha256) を作って比較し、
  不一致なら拒否する (別順序・別種集合・別 DB の場を黙って index で貼ると組成や T が壊れる)。**どちらかの署名が解決できない
  (solverConfig/species_db が無い・種名が DB に無い) ときも既定でエラー**。種を変える restart は `tools/convert_species_field.py`
  (擬似種の展開・名前で移す) を使う。`--force-species` で照合を無視できる (自己責任)。

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


def check_species_signatures(src_h5, dst_h5, force):
    """SRC/DST の隣の run 設定から化学種署名を作って照合する。不一致・解決不能は拒否 (force で警告に降格)。
    照合できたときは SRC 署名を返す (必須データセットの存在検査に使う)。"""
    from forge_species import species_signature, compare_signatures
    sig = {}
    for tag, h5 in (("SRC", src_h5), ("DST", dst_h5)):
        d = os.path.dirname(os.path.abspath(h5))
        try:
            sig[tag] = species_signature(d)
        except Exception as e:   # noqa: BLE001
            msg = f"{tag} の化学種署名が解決できない ({d}): {e}"
            if not force:
                raise SystemExit("[interp_field] REFUSED: " + msg + "  (solverConfig.yaml / species_db.yaml を隣に置くか --force-species)")
            print("[interp_field] WARNING (--force-species): " + msg)
            return None
    bad = compare_signatures(sig["SRC"], sig["DST"])
    if bad:
        msg = ("化学種署名が違う: " + "; ".join(bad) + ". 種の順序/集合/DB が違う場は index コピーできない。"
               " tools/convert_species_field.py SRC_res.h5 DST_input.h5 --meta DST/species_meta.yaml で名前により移す"
               " (どうしても index で貼るなら --force-species)。")
        if not force:
            raise SystemExit("[interp_field] REFUSED: " + msg)
        print("[interp_field] WARNING (--force-species): " + msg)
    else:
        print(f"[interp_field] species signature OK: {sig['SRC']['names'] or 'CPG (no species)'}")
    return sig["SRC"]


def check_required_datasets(h5path, sig):
    """署名が要求する保存量が SRC にあるか (res は原始量で代替可: Y{s}, Xi, k/omega)。無ければ拒否。"""
    from forge_species import required_conserved
    if sig is None:
        return
    with h5py.File(h5path, "r") as f:
        V = f["VALUE"]
        is_res = "P" in V and "Ux" in V
        missing = []
        for name in required_conserved(sig):
            if name in V:
                continue
            alt = {"roUx": "Ux", "roUy": "Uy", "roUz": "Uz", "roXi": "Xi"}.get(name)
            if name.startswith("roY"):
                alt = name[2:]
            if is_res and alt is not None and alt in V:
                continue
            if is_res and name == "roe":
                continue   # 旧 res は CPG 式で再構成 (下の警告付き経路)
            missing.append(name)
    if missing:
        raise SystemExit(f"[interp_field] REFUSED: SRC {h5path} に必須の保存量が無い: {missing} (種数 {len(sig['names'])}, tracer {sig['tracer']})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--force-species", action="store_true",
                    help="SRC/DST の physProp.species が違っても index で貼る (通常は convert_species_field.py を使う)")
    a = ap.parse_args(); g = a.gamma

    # 化学種署名の照合 (名前・順序・MW・NASA-9 係数・温度区切り・datum・tracer)。解決不能も既定でエラー (codex 2026-09-16 M3 / result-2 M2)。
    src_sig = check_species_signatures(a.src, a.dst, a.force_species)

    check_required_datasets(a.src, src_sig)

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
            if "gammaTr" in V and "reTheta" in V:   # 遷移モデル (LM2009): 原始量 → 保存量
                fields["roGamma"] = ro*np.array(V["gammaTr"]); fields["roReth"] = ro*np.array(V["reTheta"])
            for key in V:                      # scalar transport Y* -> roY*
                if key.startswith("Y") and key[1:].isdigit():
                    fields["ro"+key] = ro*np.array(V[key])
                if key == "roXi":              # 受動トレーサ (physProp.tracer)
                    fields["roXi"] = np.array(V[key])
                # 凝縮モーメント (原始 g_<s>,Q0_<s>.. → 保存 rog_<s>,roQ0_<s>..)。2026-08-17
                if key.startswith(("g_", "Q0_", "Q1_", "Q2_")):
                    fields["ro"+key] = ro*np.array(V[key])
            if "roXi" not in fields and "Xi" in V:   # 原始量 Xi しか無い res: ρ·Xi で転送 (codex result-5 M1)
                fields["roXi"] = ro*np.clip(np.array(V["Xi"], dtype=np.float64), 0.0, 1.0)
        else:                                  # input (conserved)
            fields = {n: np.array(V[n]) for n in
                      ["ro","roUx","roUy","roUz","roe","roK","roOmega","roGamma","roReth"] if n in V}
            for key in V:
                if key.startswith("roY") and key[3:].isdigit():
                    fields[key] = np.array(V[key])
                if key.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or key == "roXi":
                    fields[key] = np.array(V[key])

    # 検査で要求した保存量が転送配列に揃っているか (書込み前; codex result-5 M1)
    if src_sig is not None:
        from forge_species import required_conserved
        lack = [n for n in required_conserved(src_sig) if n not in fields]
        if lack:
            raise SystemExit(f"[interp_field] REFUSED: 転送配列に必須の保存量が無い: {lack}")

    tree = cKDTree(cs)
    with h5py.File(a.dst, "r+") as d:
        cd = centroids(d)
        _, idx = tree.query(cd)
        moved = []
        for name, arr in fields.items():
            ds = "VALUE/"+name
            if ds in d and name != "wall_dist":
                d[ds][...] = arr[idx].astype(d[ds].dtype); moved.append(name)
            elif name.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or (name.startswith("roY") and name[3:].isdigit()) or name in ("roXi", "roGamma", "roReth"):
                # 凝縮モーメントと化学種は convert 直後の入力 h5 に無いので新規作成する (無ければ forge は第 1 種以外を 0 に
                # 初期化し、carrier では rog<=roY_w のクランプで液相が消える: codex 2026-09-16 result M2)。
                # forge は VALUE/<consName> が存在すれば読む (無ければ 0 = dry restart)。2026-08-18 / 2026-09-16
                d.create_dataset(ds, data=arr[idx].astype(d["VALUE/ro"].dtype)); moved.append(name+"(new)")
        print(f"interp {a.src} -> {a.dst}: {len(cd)} dst cells, moved {moved} (wall_dist kept)")


if __name__ == "__main__":
    main()
