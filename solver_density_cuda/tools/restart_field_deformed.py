#!/usr/bin/env python3
"""節点番号・接続が同じで**一部の節点の座標だけが違う**格子への restart (局所変形した格子 B の初期場)。

    python3 solver_density_cuda/tools/restart_field_deformed.py SRC_res.h5 DST_input.h5 \\
        [--src-mesh SRC_input.h5] [--list-moved moved.csv] [--keep-src-dtype] [--dry-run] \\
        [--gamma 1.4] [--force-species] [--forge BIN] [--dst-run DIR]

用途 (plan convection-zero-thickness-edge-reconstruction §4.2「初期場」、codex diagnose 2026-10-08 te-grid-kink):
SERN 3D の格子 B (`mesh3d.te_wake_blend_H`: 後縁下流の中間線だけを局所変形、x・z 配列・接続・節点数は A と同じ) に、
格子 A の場 (例 `run_1055_r7b_m10_A_c` の最終 res) を移す。

  * **座標が変わらない節点**: 保存量を **index コピー** (restart_field.py と同じ。写した後に SRC とビット一致を検査)。
  * **座標が動いた節点だけ**: interp_field.py と同じ方式 (SRC の最寄り節点の**原始量**から保存量を組み直す:
    ρ·U、roe は res の保存量、ρ·k・ρ·ω、ρ·Y、roXi、凝縮モーメント、遷移) で埋める。
    移動した節点に元の番号の値を index コピーしない (位置が違う = 同一初期場ではない。codex 2026-10-08)。
  * `wall_dist` は**常に DST のもの** (新しい格子の壁距離)。`/AUX` (処置の重み w など) は触らない = DST のもの
    (格子が違うので、使うなら DST の格子で作り直したものを使う。plan §4 の前処理の契約)。

**使い分け**:
  - 同一格子 (座標がビット一致) の継続 → `restart_field.py` (本ツールでも結果は同じだが、移した量の集合が違いうる)。
  - 番号・接続が同じで座標の一部だけ違う (局所変形) → **本ツール**。
  - 節点数・接続が違う (解像度・トポロジの変更) → `interp_field.py` (cross-mesh、全節点を最近傍補間)。本ツールは拒否する。
  - x の station を写しただけの格子 (`L_sw_exact`) で、動いた節点にも番号で写したいときは `case/46.sern_design/r7b_index_restart.py`。

検査 (書く前):
  1. 節点数が SRC・DST で同じ、ヘキサ (可視化 primal) の接続が同一 (SRC = res の MESH/CONNE か入力の VIZMESH/CONNE、
     DST = VIZMESH/CONNE)。`--src-mesh` (SRC の入力 h5) を渡すと、その座標が SRC と一致し、境界面の接続 (タグ別) が DST と
     同一であることも見る (双子節点の番号の取り違えの検出)。
  2. 動いた節点の最寄りの SRC 節点が**座標一致の双子** (カウル上流のスリットなど) なら拒否する (どちら側の値かが不定)。
  3. 化学種: 属性 (`species_hash` ほか) の扱いは restart_field と同じ (`forge_species.plan_inherit` で書く前に継承を決め、
     書いた後に `commit_inherit`。未検証の SRC と TP の宛先は既定で停止、許可は `FORGE_ALLOW_UNVERIFIED_SPECIES=1` か
     `--force-species`)。加えて interp_field の署名照合 (隣の設定から種名・順序・MW) と必須保存量の検査を通す。ただし署名が
     解決できないとき (SERN の lump は構成種の MW が内蔵表に無く、interp_field.py も同じ理由で止まる) は、plan_inherit が
     互換性ハッシュの一致を確かめていれば通し、そうでなければ上の許可があるときだけ警告して通す。

出力: 動いた節点の数・座標範囲・最大移動量 (成分別)・最寄り距離の最大、量ごとの扱い、`VERDICT: OK`。
"""
import argparse
import os
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forge_species as fsp  # noqa: E402
from interp_field import centroids, check_required_datasets, is_3d  # noqa: E402

TOOL = "restart_field_deformed"
KEEP_FROM_DST = {"wall_dist"}


def interp_fields(V, g):
    """interp_field.py の main() (SRC が res か入力 h5) と**同じ式**で、SRC の全節点の転送用保存量を作る。
    戻り = {保存量名: SRC 長の配列}。式を変えるときは interp_field.py と一緒に変えること (試験で突き合わせている)。"""
    if "P" in V and "Ux" in V:            # res (primitives)
        ro = np.array(V["ro"]); P = np.array(V["P"])
        Ux = np.array(V["Ux"]); Uy = np.array(V["Uy"]); Uz = np.array(V["Uz"])
        if "roe" in V:
            roe = np.array(V["roe"])
        else:
            print(f"[{TOOL}] WARNING: SRC に roe が無い -> CPG 式 P/(γ-1)+½ρu² (γ = {g}) で組み直す")
            roe = P / (g - 1.0) + 0.5 * ro * (Ux**2 + Uy**2 + Uz**2)
        fields = {"ro": ro, "roUx": ro * Ux, "roUy": ro * Uy, "roUz": ro * Uz, "roe": roe}
        if "k" in V and "omega" in V:
            fields["roK"] = ro * np.array(V["k"]); fields["roOmega"] = ro * np.array(V["omega"])
        if "roGamma" in V and "roReth" in V:
            fields["roGamma"] = np.array(V["roGamma"]); fields["roReth"] = np.array(V["roReth"])
        elif "gammaTr" in V and "reTheta" in V:
            fields["roGamma"] = ro * np.array(V["gammaTr"]); fields["roReth"] = ro * np.array(V["reTheta"])
        for key in V:
            if key.startswith("Y") and key[1:].isdigit():
                fields["ro" + key] = ro * np.array(V[key])
            if key == "roXi":
                fields["roXi"] = np.array(V[key])
            if key.startswith(("g_", "Q0_", "Q1_", "Q2_")):
                fields["ro" + key] = ro * np.array(V[key])
        if "roXi" not in fields and "Xi" in V:
            fields["roXi"] = ro * np.clip(np.array(V["Xi"], dtype=np.float64), 0.0, 1.0)
    else:                                  # input (conserved)
        fields = {n: np.array(V[n]) for n in
                  ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roGamma", "roReth"] if n in V}
        for key in V:
            if key.startswith("roY") and key[3:].isdigit():
                fields[key] = np.array(V[key])
            if key.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or key == "roXi":
                fields[key] = np.array(V[key])
    return fields


def _creatable(name):
    """DST に無ければ作ってよい量 (interp_field と同じ: 化学種・凝縮モーメント・トレーサ・遷移)。"""
    return (name.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_")) or (name.startswith("roY") and name[3:].isdigit())
            or name in ("roXi", "roGamma", "roReth"))


def species_signatures(a, species_plan):
    """interp_field の署名照合 (隣の solverConfig.yaml・species_db.yaml から種名・順序・MW ほか)。戻り = SRC の署名か None。

    - plan_inherit が SRC の記録と宛先の互換性ハッシュ (`forge --resolve-species`) の一致を確かめた (species_plan あり) なら、
      署名が解決できなくても通す (ハッシュが種名・順序・係数を含む。restart_field と同じ扱い)。
    - 解決できて違えば拒否 (`--force-species` で警告に降格)。
    - 解決できない (SERN の lump の構成種は内蔵表に MW が無く、interp_field.py も同じ理由で止まる) ときは、
      `--force-species` か `FORGE_ALLOW_UNVERIFIED_SPECIES=1` (その実行だけの許可) なら警告して通し、それ以外は拒否する。"""
    from forge_species import compare_signatures, species_signature
    sig = {}
    for tag, h5 in (("SRC", a.src), ("DST", a.dst)):
        dr = a.dst_run if (tag == "DST" and a.dst_run) else os.path.dirname(os.path.abspath(h5))
        try:
            sig[tag] = species_signature(dr)
        except Exception as e:  # noqa: BLE001
            msg = f"{tag} の化学種署名が解決できない ({dr}): {e}"
            if species_plan is not None:
                print(f"[{TOOL}] {msg} — 互換性ハッシュは plan_inherit で一致を確認済みなので続ける")
                return None
            if a.force_species or fsp.allow_unverified_species():
                print(f"[{TOOL}] WARNING ({'--force-species' if a.force_species else 'FORGE_ALLOW_UNVERIFIED_SPECIES=1'}): {msg}")
                return None
            sys.exit(f"[{TOOL}] REFUSED: {msg}  (隣に solverConfig.yaml / species_db.yaml を置くか、同じ種構成と確かめたうえで"
                     " FORGE_ALLOW_UNVERIFIED_SPECIES=1 か --force-species)")
    diffs = compare_signatures(sig["SRC"], sig["DST"])
    unv = [x for x in diffs if "unverifiable" in x]       # 内蔵種の係数が設定から分からない「照合不能」(interp_field と同じ扱い)
    bad = [x for x in diffs if "unverifiable" not in x]
    if unv and not bad and not a.force_species and species_plan is None:
        if not fsp.allow_unverified_species():
            sys.exit(f"[{TOOL}] REFUSED (nothing written): " + "; ".join(unv) + " — SRC is unverified.\n" + fsp.UNVERIFIED_GUIDANCE)
        print(f"[{TOOL}] WARNING: " + "; ".join(unv) + " — allowed for this invocation by FORGE_ALLOW_UNVERIFIED_SPECIES=1")
    if bad:
        msg = "化学種署名が違う: " + "; ".join(bad)
        if not a.force_species:
            sys.exit(f"[{TOOL}] REFUSED: {msg} (種の順序/集合/DB が違う場は番号で移せない。tools/convert_species_field.py)")
        print(f"[{TOOL}] WARNING (--force-species): {msg}")
        return None
    print(f"[{TOOL}] species signature OK: {sig['SRC']['names'] or 'CPG (no species)'}")
    return sig["SRC"]


def hex_conne(f):
    """primal ヘキサの接続 (XDMF 型 9 の列)。入力 h5 は VIZMESH/CONNE、res は MESH/CONNE。無ければ None。"""
    for key in ("VIZMESH/CONNE", "MESH/CONNE"):
        if key in f:
            c = np.asarray(f[key][...], dtype=np.int64)
            if c.size and c.size % 9 == 0 and np.all(c[::9] == 9):
                return c
    return None


def bface_sets(f):
    if "BCONDS" not in f:
        return None
    out = {}
    for k in f["BCONDS"]:
        g = f["BCONDS/" + k]
        if "vizBfaceSizes" in g and "vizBfaceNodes" in g:
            out[int(k)] = (np.asarray(g["vizBfaceSizes"][...]), np.asarray(g["vizBfaceNodes"][...]))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src", help="継続元の res_*.h5 (格子 A)")
    ap.add_argument("dst", help="格子 B の forge 入力 h5 (prepare が作った sern.h5 等)。/VALUE を上書きする")
    ap.add_argument("--src-mesh", help="SRC の入力 h5 (格子 A の sern.h5)。座標・境界面の接続の照合に使う")
    ap.add_argument("--list-moved", help="動いた節点の一覧 CSV (id, DST x,y,z, Δx,Δy,Δz, 最寄り SRC id, 距離) を書く")
    ap.add_argument("--keep-src-dtype", action="store_true",
                    help="DST のデータセットを SRC の型で作り直す (restart_field.py と同じ。倍精度 res → 倍精度の種)")
    ap.add_argument("--dry-run", action="store_true", help="検査と集計だけ (書かない)")
    ap.add_argument("--gamma", type=float, default=1.4, help="SRC に roe が無い旧 res のときだけ使う (interp_field と同じ)")
    ap.add_argument("--force-species", action="store_true", help="化学種の照合を外す (属性は付けない = 未検証のまま)")
    ap.add_argument("--forge", help="--resolve-species を持つ forge (既定: FORGE_BIN, solver_density_cuda/build/forge)")
    ap.add_argument("--dst-run", help="宛先 run ディレクトリ (solverConfig.yaml の場所; 既定: DST h5 の隣)")
    a = ap.parse_args(argv)

    # --- 化学種 (書き込み前)。属性の扱いは restart_field と同じ (plan_inherit → 書いた後に commit_inherit) ---
    try:
        species_plan = fsp.plan_inherit(a.src, a.dst_run or os.path.dirname(os.path.abspath(a.dst)), forge=a.forge,
                                        force=a.force_species, tool=TOOL, inplace=not a.dry_run)
    except fsp.SpeciesCheckError as e:
        sys.exit(f"[{TOOL}] REFUSED (nothing written): {e}")
    src_sig = species_signatures(a, species_plan)
    check_required_datasets(a.src, src_sig)

    # --- 格子の照合 ---
    with h5py.File(a.src, "r") as s, h5py.File(a.dst, "r") as d:
        if "VALUE" not in s or "VALUE" not in d:
            sys.exit(f"[{TOOL}] REFUSED: SRC / DST のどちらかに /VALUE が無い")
        n_src, n_dst = s["VALUE/ro"].shape[0], d["VALUE/ro"].shape[0]
        cs_all = np.asarray(s["MESH/COORD"][...], dtype=np.float64).reshape(-1, 3)
        cd_all = np.asarray(d["MESH/COORD"][...], dtype=np.float64).reshape(-1, 3)
        if not (n_src == n_dst == cs_all.shape[0] == cd_all.shape[0]):
            sys.exit(f"[{TOOL}] REFUSED: 節点数が違う (SRC 値 {n_src} / 座標 {cs_all.shape[0]}, DST 値 {n_dst} / 座標 {cd_all.shape[0]})。"
                     " 番号の対応が無い格子は interp_field.py を使うこと")
        hs, hd = hex_conne(s), hex_conne(d)
        if hs is None or hd is None:
            sys.exit(f"[{TOOL}] REFUSED: ヘキサの接続 (MESH/CONNE / VIZMESH/CONNE) が {'SRC' if hs is None else 'DST'} に無い (接続の同一性を確かめられない)")
        if not np.array_equal(hs, hd):
            sys.exit(f"[{TOOL}] REFUSED: ヘキサの接続が SRC と DST で違う。番号の対応が無いので interp_field.py を使うこと")
        for tag, perm_f in (("SRC", s), ("DST", d)):
            if "MESH/RENUMBER_PERM" in perm_f:
                print(f"[{TOOL}] {tag} は節点を再番号付け済み (MESH/RENUMBER_PERM)。接続が一致したので同じ番号として扱う")
        n_bf = None
        if a.src_mesh:
            with h5py.File(a.src_mesh, "r") as m:
                cm = np.asarray(m["MESH/COORD"][...], dtype=np.float64).reshape(-1, 3)
                if cm.shape != cs_all.shape or not np.array_equal(cm, cs_all):
                    sys.exit(f"[{TOOL}] REFUSED: --src-mesh の座標が SRC (res) の MESH/COORD と一致しない (別の格子の res)")
                hm = hex_conne(m)
                if hm is not None and not np.array_equal(hm, hd):
                    sys.exit(f"[{TOOL}] REFUSED: --src-mesh と DST のヘキサの接続が違う")
                bs, bd = bface_sets(m), bface_sets(d)
                if bs is None or bd is None or set(bs) != set(bd) or not all(
                        np.array_equal(bs[k][0], bd[k][0]) and np.array_equal(bs[k][1], bd[k][1]) for k in bs):
                    sys.exit(f"[{TOOL}] REFUSED: 境界面の接続 (BCONDS/*/vizBface*) が --src-mesh と DST で違う (双子の取り違えの恐れ)")
                n_bf = len(bs)

    diff = cd_all - cs_all
    moved = np.where(np.any(diff != 0.0, axis=1))[0]
    still = np.where(~np.any(diff != 0.0, axis=1))[0]
    print(f"[{TOOL}] 節点 {n_dst}: 接続同一" + (f"・境界面の接続同一 ({n_bf} 群)" if n_bf is not None else "")
          + f"、座標が一致 {still.size}、**動いた {moved.size}**")
    if moved.size == 0:
        print(f"[{TOOL}] 座標が全節点で一致 = 同一格子。restart_field.py を使うのが本来 (このまま index コピーで続ける)")
    else:
        lo, hi = cd_all[moved].min(0), cd_all[moved].max(0)
        mx = np.abs(diff[moved]).max(0)
        print(f"[{TOOL}] 動いた節点の範囲 (DST): x {lo[0]:.6g}..{hi[0]:.6g}  y {lo[1]:.6g}..{hi[1]:.6g}  z {lo[2]:.6g}..{hi[2]:.6g}")
        print(f"[{TOOL}] 最大移動量: |Δx| {mx[0]:.4e}  |Δy| {mx[1]:.4e}  |Δz| {mx[2]:.4e}")

    # --- 転送する量 (interp_field と同じ集合) ---
    with h5py.File(a.src, "r") as s:
        V = s["VALUE"]
        fields = interp_fields(V, a.gamma)
        src_cons = {n: np.array(V[n]) for n in fields if n in V}
        s3 = is_3d(s)
        with h5py.File(a.dst, "r") as d0:
            d3 = is_3d(d0)
            dst_names = set(d0["VALUE"].keys())
        nd = 3 if (s3 and d3) else 2
        cs = centroids(s, nd)
    if src_sig is not None:
        from forge_species import required_conserved
        lack = [n for n in required_conserved(src_sig) if n not in fields]
        if lack:
            sys.exit(f"[{TOOL}] REFUSED: 転送配列に必須の保存量が無い: {lack}")
    # 不動節点は SRC の保存量そのものを写す。保存量の無い量は index コピーできないので止める
    # (restart_field は「SRC に無く据置」で DST の初期値を残すが、ここでは動いた節点との整合が崩れるので通さない)
    no_cons = [n for n in fields if n not in src_cons and (n in dst_names or _creatable(n))]
    if no_cons and still.size:
        sys.exit(f"[{TOOL}] REFUSED: SRC に保存量が無く不動節点を index コピーできない量: {no_cons}")
    extra = sorted(n for n in dst_names - KEEP_FROM_DST if n not in fields)
    with h5py.File(a.src, "r") as s:
        extra_in_src = [n for n in extra if n in s["VALUE"]]
    if extra_in_src:
        sys.exit(f"[{TOOL}] REFUSED: DST にあり SRC にもあるが補間の式が無い量 {extra_in_src} (動いた節点を埋められない)")

    # --- 動いた節点の最寄り (interp_field と同じ木・同じ照合座標) ---
    idx = dist = None
    if moved.size:
        tree = cKDTree(cs)
        with h5py.File(a.dst, "r") as d0:
            cd = centroids(d0, nd)
        dist, idx = tree.query(cd[moved])
        # 最寄りの SRC 節点が座標一致の双子 (自分以外に距離 0 の節点がある) なら、どちら側の値かが不定
        dd2, _ = tree.query(cs[idx], k=2)
        twin = dd2[:, 1] == 0.0
        if np.any(twin):
            bad = moved[twin][:5]
            sys.exit(f"[{TOOL}] REFUSED: 動いた節点 {int(twin.sum())} 個の最寄りの SRC 節点が座標一致の双子 (例: DST {bad.tolist()})。"
                     " 壁の両側を取り違える恐れがある")
        k = int(np.argmax(dist))
        print(f"[{TOOL}] 最寄りの SRC 節点までの距離: max {float(dist.max()):.4e} (DST {int(moved[k])} -> SRC {int(idx[k])})、"
              f"median {float(np.median(dist)):.4e}、距離 0 {int((dist == 0).sum())}、使った SRC 節点 {np.unique(idx).size}、双子 0")
        if a.list_moved:
            np.savetxt(a.list_moved, np.column_stack([moved, cd_all[moved], diff[moved], idx, dist]), delimiter=",",
                       header="dst_id,x,y,z,dx,dy,dz,src_id,dist", comments="",
                       fmt=["%d", "%.9g", "%.9g", "%.9g", "%.6e", "%.6e", "%.6e", "%d", "%.6e"])
            print(f"[{TOOL}] 動いた節点の一覧: {a.list_moved}")

    names = [n for n in fields if (n in dst_names or _creatable(n)) and n not in KEEP_FROM_DST]
    created = [n for n in names if n not in dst_names]
    print(f"[{TOOL}] 移す量: {names}" + (f"  (DST に無く新規作成: {created})" if created else "")
          + f"  /  DST を残す: {sorted(KEEP_FROM_DST & dst_names)}")
    if a.dry_run:
        print(f"[{TOOL}] --dry-run: 何も書いていない")
        return 0

    fsp.write_species_attrs(a.dst, None)      # 書き込み途中で失敗しても古い属性が残らないように先に消す
    narrowed, bad = [], []
    with h5py.File(a.dst, "r+") as d:
        for n in names:
            ds = "VALUE/" + n
            src_arr = src_cons.get(n)
            if ds in d:
                if a.keep_src_dtype and src_arr is not None and d[ds].dtype != src_arr.dtype:
                    old = np.asarray(d[ds][...]); del d[ds]
                    d.create_dataset(ds, data=old.astype(src_arr.dtype))
            else:
                d.create_dataset(ds, data=np.zeros(n_dst, dtype=(src_arr.dtype if (a.keep_src_dtype and src_arr is not None)
                                                                   else d["VALUE/ro"].dtype)))
            dt = d[ds].dtype
            out = np.asarray(d[ds][...])
            if still.size:
                out[still] = src_arr[still].astype(dt, copy=False)
            if moved.size:
                out[moved] = fields[n][idx].astype(dt)
            d[ds][...] = out
            # 検査: 不動節点は SRC とビット一致 (型の縮小は大きさを報告)、動いた節点は補間値と一致
            got = np.asarray(d[ds][...])
            if still.size and not np.array_equal(got[still], src_arr[still]):
                if src_arr.dtype.itemsize > got.dtype.itemsize:
                    rel = np.max(np.abs(got[still].astype(np.float64) - src_arr[still].astype(np.float64))) / max(np.max(np.abs(src_arr)), 1e-300)
                    narrowed.append((n, rel))
                else:
                    bad.append(n + " (不動節点)")
            if moved.size and not np.array_equal(got[moved], fields[n][idx].astype(dt)):
                bad.append(n + " (動いた節点)")
        if bad:
            sys.exit(f"[{TOOL}] 検査 NG: 書いた値が一致しない: {bad}")
        fsp.commit_inherit(d, species_plan)
    if narrowed:
        print(f"[{TOOL}] ⚠ 型の縮小で丸めが入った (倍精度の場を種にするなら --keep-src-dtype):")
        for n, rel in narrowed:
            print(f"    {n:<10} 相対 {rel:.3e}")
    print(f"[{TOOL}] species 属性: {('継承 (species_input_unverified=%d)' % species_plan['species_input_unverified']) if species_plan else 'なし (未検証のまま)'}")
    print(f"VERDICT: OK ({len(names)} 量: 不動節点 {still.size} は index コピー"
          + ("で SRC とビット一致" if not narrowed else f" (うち {len(narrowed)} 量は型の縮小で丸めあり)")
          + f"、動いた節点 {moved.size} は interp_field と同じ原始量の最近傍補間; wall_dist・/AUX は DST のまま)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
