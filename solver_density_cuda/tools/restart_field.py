#!/usr/bin/env python3
"""同一メッシュの restart: res_*.h5 の**保存量をそのまま**入力 h5 の /VALUE へ写す。

    python3 solver_density_cuda/tools/restart_field.py SRC_res.h5 DST_input.h5 [--dry-run] [--src-mesh SRC_MESH.h5]

**なぜ専用ツールが要るか** (2026-09-23 に実際にやらかした):
`res_*.h5` は保存量 (`ro,roUx,roUy,roUz,roe,roK,roOmega,…`) と原始量 (`P,T,Ux,…`) の**両方**を持つ。
同一メッシュの restart で原始量から保存量を**組み直すと**、`roUx = ro * Ux` の丸めで元の値に戻らない。
`case/56.gap_tp1187` で `interp_field.py` (cross-mesh 用) を同一メッシュに使ったところ、
`roUy` が **62104/65194 セルで不一致・最大 1.6 %** になった。測る量が
$\\dot m = |\\int \\rho U_y dx|$ だったので、そのまま回せば実験ごと壊れていた
(既存の `case/*/restart_field.py` も原始量から組み直す実装で、同じ罠を持つ)。

本ツールは**保存量を index コピー**する。DST/VALUE にある各データセットについて、
SRC/VALUE に同名があればその値をそのまま書き、無ければ DST の値を残す。
`wall_dist` は**常に DST のものを残す** (メッシュ由来の量)。
最後に「SRC と DST が保存量でビット一致すること」を検査し、しなければ**失敗させる**。

メッシュが違うとき (解像度変更・quad↔tri) は本ツールでなく `interp_field.py` を使う。

**同一メッシュの検査** (plan tooling-generator-defaults-and-launch-gate §4.5・§5.1 #20、2026-10-15): 書き込みの前に、
SRC と DST がそれぞれ載っている**格子** (変換済みの forge 入力 h5) を決め、両者が同一であることを確かめる。確かめられなければ
**何も書かずに停止**する。2026-10-15 まで接続を見ず、座標が無くても続行し、倍精度の座標を float32 に丸めて比べていたので、
座標を約 1e−8 変えた別の入力も「同一メッシュ」と判定した (codex diagnose `notes/reviews/2026-10-15-gate-profile-judgments-diagnose.md`)。

1. **格子を決める**: 変換済みの h5 (`/BCONDS`・`/PLANES`・`/CELLS`・`/DUAL`・`/VIZMESH` のどれかを持つ) はそれ自身が格子。
   `res_*.h5` (forge の出力は `/MESH/COORD`・`/MESH/CONNE`・`/VALUE` だけ: `output.cpp:149`・`:191`) は格子を持たないので、
   SRC は `--src-mesh` か、SRC の隣の `solverConfig.yaml` の `mesh.meshFileName`、DST は `--dst-run` (既定 DST の隣) の
   `solverConfig.yaml` の `mesh.meshFileName` を格子とする。決められない (config が無い・読めない・書いていない・無いファイル・
   格子の構造を持たない) ときは停止する。
2. **格子を持たない側を格子に結び付ける**: その h5 の `/MESH/COORD`・`/MESH/CONNE` が必須。座標は**その h5 の格納の型**で
   格子の座標と完全一致 (float のビルドの res は float32 で書くので、倍精度の格子は float32 に丸めて比べる。res はそれ以上の桁を
   持たない)、接続は node なら格子の `/VIZMESH/CONNE`、cell なら `/MESH/CONNE` と値が完全一致 (forge の出力が書く接続)。
3. **SRC の格子と DST の格子が同一** (同じファイルなら省略): `mesh_provenance.mesh_parts` の部分 `coord` (座標を**格納の型のまま**
   バイト単位で。型が違えば別の入力)・`topology` (接続・CV の順・再番号・個数)・`boundary` (境界の所属と `bcondKind`)・
   `geometry` (面ベクトル・体積・重心)・`other` (分類していない位置) のハッシュがすべて一致。`wall_dist` は DST のものを残すので
   比べず、違えば注記だけ出す。どちらの格子も `/MESH/COORD`・`/MESH/CONNE` が無ければ停止する。
4. **座標の精度だけを上げた再変換への引き継ぎ** (`--allow-coord-upcast`、plan §5.1 #20b・#20a (4)、2026-10-15): 3 は座標の型だけ違う
   同じ `.msh` の再変換 (float32 の変換器 → 倍精度の変換器) を別の格子として拒む。旧い float32 の格子の場を倍精度の格子へ移すときだけ、
   明示してこの逃げ道を使う。条件: `topology`・`boundary`・`other` が同一、座標は SRC の格子の方が狭い浮動小数の型で、DST の格子の座標を
   その型に丸めると SRC の格子の座標と**完全一致**、`geometry` は違ってよい (同じ座標から計算した精度の差; 最大の相対差を報告)。
   それ以外 (型が同じ・広げる向きが逆・丸めても違う座標・接続や境界の違い) は 3 と同じく停止する。通ったときは VERDICT の行に
   「座標の精度の引き上げ (旧 → 新 = 格子の不変部分の sha256、`mesh_provenance` の mesh-parts)」を書き、DST の h5 のルートの属性
   `restart_coord_upcast` (JSON。`mesh_provenance.COORD_UPCAST_ATTR`; 格子の不変部分のハッシュには入らない) に記録する。投入の関門
   (`launch_gate.py`) がその記録を `LAUNCH_GATE.txt` と `RUN_PROVENANCE.txt` に写す (`run_case.sh` は起動ごとに `RUN_PROVENANCE.txt` を
   書き直すので、ここから直接は書かない)。旧い場の保存量は新しい双対の体積で読み直されるので、引き継いだ run は段の判定で
   「restart の跳ね」を見る (skill forge-aws-run §3)。この記録は旗なしの restart (同じ格子の場で上書き) で消す。

**化学種の属性** (plans/active/thermophysics-solver-owned-species-db.md §4.3, #3b): 書き込み前に SRC の属性
(`species_hash` ほか) と SRC の隣の解決済み記録 (完全性ハッシュ再計算) を検証し、宛先 run (`--dst-run`, 既定は DST の隣) を
`forge --resolve-species` (`--forge` / `FORGE_BIN`) で解決して互換性ハッシュが一致したときだけ属性を DST に継承する
(記録も DST の隣へ複製)。不一致は差のある係数を示して**書き込まずに停止** (`--force-species` で属性なしのまま通す)。
SRC が未検証 (属性なし / `species_input_unverified=1`) で宛先が TP (または solverConfig.yaml が無く判定できない) とき・
宛先を解決できない (旧バイナリ) ときは**既定で書き込まずに停止** (ソルバと同じ規約, #3c)。許可はその実行だけの
`FORGE_ALLOW_UNVERIFIED_SPECIES=1` か `--force-species` で、そのとき DST の属性は消す (宛先のハッシュで埋めない;
ソルバ側でも未検証として扱われ、その run にも同じ許可が要る)。CPG は対象外。
ただし印付きの SRC (`species_input_unverified=1`) で `species_hash` = 宛先ハッシュなら、ソルバ (入力場照合の一致分岐) と同じく
許可なしで通し、DST に同じハッシュと印を継承する (#3d; 印は消さない)。ハッシュ不一致の印付きは上と同じく停止。
"""
import argparse, os, sys
import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forge_species as fsp  # noqa: E402
import mesh_provenance as mpv  # noqa: E402

KEEP_FROM_DST = {"wall_dist"}
# 変換済みの格子だけが持つグループ (forge の res は /MESH/COORD・/MESH/CONNE・/VALUE・/CHECKPOINT だけを書く: output.cpp)
MESH_GROUPS = ("BCONDS", "PLANES", "CELLS", "DUAL", "VIZMESH")
# mesh_parts の部分のうち SRC と DST で同一を求めるもの (wall_dist は DST のものを残すので比べない)
SAME_PARTS = ("coord", "topology", "geometry", "boundary", "other")


class MeshIdentityError(Exception):
    """同一メッシュと確かめられない (格子を決められない・座標や接続が無い・一致しない)。"""


def _carries_mesh(f):
    return any(g in f for g in MESH_GROUPS)


def _same_file(x, y):
    try:
        return os.path.samefile(x, y)
    except OSError:
        return os.path.abspath(x) == os.path.abspath(y)


def _run_mesh_file(run_dir, what):
    """run の solverConfig.yaml の mesh.meshFileName → 格子のパス (forge と同じく run ディレクトリからの相対)。"""
    import yaml_strict
    cfg_path = os.path.join(run_dir, "solverConfig.yaml")
    if not os.path.isfile(cfg_path):
        raise MeshIdentityError(f"{what} の格子を決められない: {cfg_path} が無い")
    try:
        with open(cfg_path, encoding="utf-8") as fh:
            doc = yaml_strict.load(fh.read())
    except Exception as e:                          # 重複キー・merge key・構文の誤り (yaml-cpp と解釈が割れうる)
        raise MeshIdentityError(f"{what} の格子を決められない: {cfg_path} を読めない ({e})") from None
    mesh = (doc or {}).get("mesh") if isinstance(doc, dict) else None
    name = mesh.get("meshFileName") if isinstance(mesh, dict) else None
    if not isinstance(name, str) or not name.strip():
        raise MeshIdentityError(f"{what} の格子を決められない: {cfg_path} に mesh.meshFileName が無い")
    return os.path.join(run_dir, name.strip()), f"{cfg_path} の mesh.meshFileName"


def _resolve_mesh(h5, what, explicit, run_dir, flag):
    """h5 が載っている格子を決める → (格子のパス, 決め方)。h5 自身が格子なら (h5, '自身')。"""
    with h5py.File(h5, "r") as f:
        own = _carries_mesh(f)
    if own:
        if explicit:
            raise MeshIdentityError(f"{what} {h5} は格子を持つので {flag} は使わない (格子を持たない res のときだけ)")
        return h5, "自身"
    if explicit:
        mesh, how = explicit, flag
    else:
        mesh, how = _run_mesh_file(run_dir, f"{what} {os.path.basename(h5)} (格子を持たない)")
        how += f" (別の格子なら {flag})" if flag else ""
    if not os.path.isfile(mesh):
        raise MeshIdentityError(f"{what} の格子 {mesh} が無い ({how})")
    with h5py.File(mesh, "r") as f:
        if not _carries_mesh(f):
            raise MeshIdentityError(f"{what} の格子 {mesh} ({how}) は変換済みの格子でない "
                                    f"({'・'.join('/' + g for g in MESH_GROUPS)} のどれも無い)")
    return mesh, how


def _need_coord_conne(f, path, what):
    miss = [k for k in ("MESH/COORD", "MESH/CONNE") if k not in f]
    if miss:
        raise MeshIdentityError(f"{what} {path} に {'・'.join('/' + k for k in miss)} が無い (座標・接続を照合できない)")


def _bind_to_mesh(h5, mesh, what):
    """格子を持たない h5 (res) の座標・接続が格子のものか。戻り = 照合の説明。"""
    with h5py.File(h5, "r") as f, h5py.File(mesh, "r") as m:
        _need_coord_conne(f, h5, what)
        _need_coord_conne(m, mesh, f"{what} の格子")
        fc, mc = f["MESH/COORD"], m["MESH/COORD"]
        if fc.dtype.kind != "f" or mc.dtype.kind != "f":
            raise MeshIdentityError(f"{what}: /MESH/COORD が浮動小数でない ({fc.dtype} / {mc.dtype})")
        if fc.shape != mc.shape:
            raise MeshIdentityError(f"{what} {os.path.basename(h5)} の座標の形 {fc.shape} が格子 {mesh} の {mc.shape} と違う")
        a, b = np.asarray(fc), np.asarray(mc)
        # h5 の格納の型で比べる: 格子の方が広い型 (float のビルドの res と倍精度の格子) なら格子をその型に丸める。
        # 格子の方が狭い型 (倍精度のビルドの res と float32 の格子) なら広げるのは厳密なのでそのまま一致を求めることになる
        if not np.array_equal(a, b.astype(a.dtype)):
            nd = int(np.count_nonzero(a != b.astype(a.dtype)))
            raise MeshIdentityError(f"{what} {os.path.basename(h5)} の座標が格子 {mesh} と {nd}/{a.size} 成分で違う "
                                    f"({a.dtype} で照合)。この格子の場でない")
        ck = "VIZMESH/CONNE" if "VIZMESH/CONNE" in m else "MESH/CONNE"   # node の出力は primal の接続 (output.cpp:154)
        fa, ma = np.asarray(f["MESH/CONNE"]), np.asarray(m[ck])
        if fa.shape != ma.shape or not np.array_equal(fa.astype(np.int64), ma.astype(np.int64)):
            raise MeshIdentityError(f"{what} {os.path.basename(h5)} の接続 /MESH/CONNE {fa.shape} が格子 {mesh} の /{ck} "
                                    f"{ma.shape} と違う。この格子の場でない")
        how = f"座標 ({a.dtype}" + (f"; 格子の {b.dtype} を {a.dtype} に丸めて" if b.dtype.itemsize > a.dtype.itemsize else "") \
              + f")・接続 (/{ck}) が一致"
    return how


def _attr_eq(x, y):
    tx, ty = (bytes(v).decode("utf-8", "replace") if isinstance(v, (bytes, np.bytes_)) else v for v in (x, y))
    if isinstance(tx, str) or isinstance(ty, str):
        return tx == ty
    return np.asarray(tx).shape == np.asarray(ty).shape and np.array_equal(np.asarray(tx), np.asarray(ty))


def _differing_members(fa, fb, ma, mb):
    """部分の中で違う位置 (案内用。判定は mesh_parts のハッシュで決まる)。"""
    sa, sb = set(ma), set(mb)
    out = [f"{k} (SRC の格子だけ)" for k in sorted(sa - sb)] + [f"{k} (DST の格子だけ)" for k in sorted(sb - sa)]
    for k in sorted(sa & sb):
        if "@" in k:
            p, an = k.split("@", 1)
            same = _attr_eq(fa[p].attrs[an], fb[p].attrs[an])
        else:
            da, db = fa[k], fb[k]
            same = da.dtype.str == db.dtype.str and da.shape == db.shape and np.array_equal(da[()], db[()])
        if not same:
            if k == "MESH/COORD":
                da, db = fa[k], fb[k]
                d = (f"{k} ({da.dtype} {da.shape} / {db.dtype} {db.shape}"
                     + (f", 最大 |Δ| {np.max(np.abs(np.asarray(da, np.float64) - np.asarray(db, np.float64))):.3e}"
                        if da.shape == db.shape and da.dtype.kind == db.dtype.kind == "f" else "") + ")")
                out.append(d)
            else:
                out.append(k)
    return out


def _coord_upcast(sm, dm, pa, pb, diff):
    """SRC の格子 sm と DST の格子 dm の違い diff (部分の名前) が「座標の精度の引き上げ」だけか (上の 4)。
    戻り (記録 dict, None) か (None, 拒む理由)。"""
    if not set(diff) <= {"coord", "geometry"}:
        return None, (f"部分 {', '.join(p for p in diff if p not in ('coord', 'geometry'))} も違う "
                      "(--allow-coord-upcast は座標の精度だけの違いに限る: topology・boundary・other は同一が要る)")
    if "coord" not in diff:
        return None, "座標が同一で幾何だけ違う (座標の精度の引き上げではない)"
    with h5py.File(sm, "r") as fa, h5py.File(dm, "r") as fb:
        ca, cb = fa["MESH/COORD"], fb["MESH/COORD"]
        if ca.dtype.kind != "f" or cb.dtype.kind != "f" or not ca.dtype.itemsize < cb.dtype.itemsize:
            return None, (f"座標の型 SRC の格子 {ca.dtype} → DST の格子 {cb.dtype} は精度の引き上げでない "
                          "(SRC の格子の方が狭い浮動小数の型のときだけ)")
        if ca.shape != cb.shape:
            return None, f"座標の形が違う ({ca.shape} / {cb.shape})"
        a, b = np.asarray(ca), np.asarray(cb)
        rb = b.astype(a.dtype)
        if not np.array_equal(a, rb):
            nd = int(np.count_nonzero(a != rb))
            return None, (f"DST の格子の座標を {a.dtype} に丸めても SRC の格子の座標と {nd}/{a.size} 成分で違う "
                          "(同じ .msh の再変換でない)")
        geo, worst = 0.0, None
        ga, gb = set(pa["members"].get("geometry", [])), set(pb["members"].get("geometry", []))
        if ga != gb:
            return None, f"幾何の位置の集合が違う ({sorted(ga ^ gb)[:4]})"
        for k in sorted(ga):
            if "@" in k:
                continue
            da, db = fa[k], fb[k]
            if da.shape != db.shape:
                return None, f"幾何 {k} の形が違う ({da.shape} / {db.shape})"
            if da.dtype.kind != "f" or db.dtype.kind != "f":
                continue
            x, y = np.asarray(da, np.float64), np.asarray(db, np.float64)
            scale = float(np.max(np.abs(y))) if y.size else 0.0
            r = float(np.max(np.abs(x - y))) / scale if scale > 0 else 0.0
            if r > geo:
                geo, worst = r, k
    rec = {"kind": "coord_upcast", "tool": "restart_field.py --allow-coord-upcast",
           "src_mesh": os.path.abspath(sm), "src_coord_dtype": str(ca.dtype), "src_invariant": pa["invariant"],
           "dst_mesh": os.path.abspath(dm), "dst_coord_dtype": str(cb.dtype), "dst_invariant": pb["invariant"],
           "parts_scheme": mpv.PARTS_SCHEME, "geometry_max_rel_diff": geo, "geometry_max_rel_diff_at": worst}
    return rec, None


def check_same_mesh(src, dst, src_mesh_arg=None, dst_run=None, allow_coord_upcast=False):
    """書き込みの前に SRC と DST が同一メッシュの上の場であることを確かめる。戻り = (報告の行, 座標の精度の引き上げの記録か None)。
    MeshIdentityError で停止。allow_coord_upcast は上の 4。"""
    lines = []
    sm, s_how = _resolve_mesh(src, "SRC", src_mesh_arg, os.path.dirname(os.path.abspath(src)), "--src-mesh")
    dm, d_how = _resolve_mesh(dst, "DST", None, dst_run or os.path.dirname(os.path.abspath(dst)), "")
    for h5, mesh, how, what in ((src, sm, s_how, "SRC"), (dst, dm, d_how, "DST")):
        if how == "自身":
            with h5py.File(mesh, "r") as f:
                _need_coord_conne(f, mesh, f"{what} の格子")
            lines.append(f"{what} の格子: {mesh} (自身)")
        else:
            lines.append(f"{what} の格子: {mesh} ({how}) — {os.path.basename(h5)} の {_bind_to_mesh(h5, mesh, what)}")
    if _same_file(sm, dm):
        lines.append("SRC と DST の格子は同じファイル")
        return lines, None
    pa, pb = mpv.mesh_parts(sm), mpv.mesh_parts(dm)
    diff = [p for p in SAME_PARTS if pa["parts"][p] != pb["parts"][p]]
    upcast = None
    if diff:
        why = None
        if allow_coord_upcast:
            upcast, why = _coord_upcast(sm, dm, pa, pb, diff)
        if upcast is None:
            with h5py.File(sm, "r") as fa, h5py.File(dm, "r") as fb:
                det = []
                for p in diff:
                    mem = _differing_members(fa, fb, pa["members"].get(p, []), pb["members"].get(p, []))
                    det.append(f"{p}: " + (", ".join(mem[:8]) + (f" ほか {len(mem) - 8}" if len(mem) > 8 else "") if mem else "(位置は同じで属性の並び等)"))
            raise MeshIdentityError(f"SRC の格子 {sm} と DST の格子 {dm} が違う (部分 {', '.join(diff)}):\n    " + "\n    ".join(det)
                                    + (f"\n  --allow-coord-upcast でも通さない: {why}" if why else "")
                                    + "\n  同一メッシュでないので interp_field.py を使う (同じ .msh を変換し直した格子でも、座標の型・値・接続が"
                                      "違えば別の入力。座標の精度だけを上げた再変換なら --allow-coord-upcast)")
        lines.append(f"座標の精度の引き上げ (--allow-coord-upcast): SRC の格子 {sm} (座標 {upcast['src_coord_dtype']}) → DST の格子 {dm} "
                     f"(座標 {upcast['dst_coord_dtype']})。topology・boundary・other は同一、DST の座標を {upcast['src_coord_dtype']} に丸めると"
                     f" SRC と完全一致、幾何の最大相対差 {upcast['geometry_max_rel_diff']:.3e}"
                     + (f" ({upcast['geometry_max_rel_diff_at']})" if upcast['geometry_max_rel_diff_at'] else "")
                     + f"。格子の不変部分の sha256 {upcast['src_invariant'][:16]} → {upcast['dst_invariant'][:16]} ({mpv.PARTS_SCHEME})")
        lines.append("注記: 旧い場の保存量は新しい双対の体積で読み直される。引き継いだ run は段の判定で restart の跳ねを見る")
    else:
        lines.append(f"SRC と DST の格子が一致: 座標 (格納の型のまま)・接続・CV の順・境界の所属・幾何 (mesh_parts {mpv.PARTS_SCHEME})")
    if pa["parts"]["wall_dist"] != pb["parts"]["wall_dist"]:
        lines.append("注記: wall_dist は SRC の格子と違う (DST のものを残す)")
    return lines, upcast


ap = argparse.ArgumentParser()
ap.add_argument("src", help="継続元の res_*.h5 (か変換済みの入力 h5)")
ap.add_argument("dst", help="forge 入力 h5 (mesh.h5 等)。/VALUE を上書きする")
ap.add_argument("--src-mesh", help="SRC が格子を持たない res のとき、その場の格子 (変換済みの入力 h5)。"
                                   "既定: SRC の隣の solverConfig.yaml の mesh.meshFileName。どちらも無ければ停止する")
ap.add_argument("--dry-run", action="store_true", help="書かずに何が起きるかだけ出す")
ap.add_argument("--keep-src-dtype", action="store_true",
                help="DST のデータセットを SRC の型で作り直す (倍精度 res → 倍精度 seed)。"
                     "**FP64 ビルドは倍精度の入力をそのまま読める** ので、倍精度の場を種にするときはこれを使う")
ap.add_argument("--dst-run", help="宛先 run ディレクトリ (solverConfig.yaml の場所; 既定: DST h5 の隣)")
ap.add_argument("--forge", help="--resolve-species を持つ forge (既定: FORGE_BIN, solver_density_cuda/build/forge)")
ap.add_argument("--force-species", action="store_true",
                help="化学種の不一致・記録の欠落でも保存量を写す (属性は付けない = 未検証のまま)")
ap.add_argument("--allow-coord-upcast", action="store_true",
                help="座標の精度だけを上げた同じ .msh の再変換 (float32 の格子 → 倍精度の格子) への引き継ぎを許す。"
                     "接続・境界・その他は同一、DST の座標を SRC の格子の型に丸めると完全一致のときだけ (docstring の 4)")
a = ap.parse_args()

# --- 同一メッシュの検査 (書き込み前、化学種の解決より先; plan tooling-generator-defaults-and-launch-gate §5.1 #20・#20b) ---
try:
    mesh_lines, upcast = check_same_mesh(a.src, a.dst, a.src_mesh, a.dst_run, allow_coord_upcast=a.allow_coord_upcast)
except (MeshIdentityError, OSError) as e:
    sys.exit(f"[restart_field] REFUSED (nothing written): 同一メッシュと確かめられない: {e}")
if upcast is not None:
    import datetime
    upcast.update({"src": os.path.abspath(a.src), "dst": os.path.abspath(a.dst),
                   "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds")})

# --- 化学種: SRC の記録を検証し、宛先を解決して継承できるか決める (書き込み前; §4.3) ---
try:
    species_plan = fsp.plan_inherit(a.src, a.dst_run or os.path.dirname(os.path.abspath(a.dst)), forge=a.forge,
                                    force=a.force_species, tool="restart_field", inplace=not a.dry_run)
except fsp.SpeciesCheckError as e:
    sys.exit(f"[restart_field] REFUSED (nothing written): {e}")
if not a.dry_run:
    fsp.write_species_attrs(a.dst, None)      # 書き込み途中で失敗しても古い属性が残らないように先に消す

with h5py.File(a.src, "r") as s, h5py.File(a.dst, "r" if a.dry_run else "r+") as d:
    if "VALUE" not in s or "VALUE" not in d:
        sys.exit("SRC / DST のどちらかに /VALUE が無い")
    sv, dv = s["VALUE"], d["VALUE"]

    # --- 同一メッシュであることを確かめる (違えば interp_field.py を使うべき) ---
    n_dst = dv[next(iter(dv))].shape[0]
    n_src = sv[next(iter(sv))].shape[0]
    if n_src != n_dst:
        sys.exit(f"セル数が違う (SRC {n_src} != DST {n_dst})。"
                 " 同一メッシュでないので interp_field.py を使うこと")
    # 格子の照合は書き込みの前に済ませた (check_same_mesh)。float のビルドの res (座標 float32) と倍精度の格子は、res を
    # その run の格子に float32 で結び付け、格子どうしは格納の型のまま比べる (plan architecture-float-state-double-geometry §6.18 の
    # 「DST を SRC の精度に丸めて照合」は格子どうしの比較には使わない; 1e−8 だけ違う別の入力を通していた)
    for ln in mesh_lines:
        print(f"[restart_field] {ln}")
    print(f"同一メッシュを確認 ({n_dst} cells)")

    moved, kept, missing = [], [], []
    for name in sorted(dv):
        if name in KEEP_FROM_DST:
            kept.append(name); continue
        if name not in sv:
            missing.append(name); continue
        val = np.asarray(sv[name])
        if val.shape != dv[name].shape:
            sys.exit(f"{name}: 形が違う (SRC {val.shape} != DST {dv[name].shape})")
        if not a.dry_run:
            if a.keep_src_dtype and dv[name].dtype != val.dtype:
                # DST のデータセットを SRC の型で作り直す (縮小丸めを起こさない)。
                del dv[name]
                dv.create_dataset(name, data=val)
            else:
                dv[name][...] = val.astype(dv[name].dtype, copy=False)
        moved.append(name)

    print(f"移した保存量  : {moved}")
    print(f"DST を残す    : {kept}")
    if missing:
        print(f"SRC に無く据置: {missing}")

    if a.dry_run:
        print("--dry-run: 何も書いていない")
        sys.exit(0)

    # --- 検査: 写した量が SRC とビット一致すること ---
    # **SRC の方が広い型のとき** (倍精度 run の res_*.h5 → float32 の入力 h5) は縮小丸めが入る。
    # ~~「forge の入力 h5 は float32 なので不可避」~~ **誤り** (2026-09-24, codex result M3):
    # 読込先は `std::vector<geom_float>` (`variables.cpp:730-736`) で **FP64 ビルドでは double**、
    # HDF5 にも制約は無い。**倍精度の場を種にするときは `--keep-src-dtype` を使うこと**。
    # これを怠って倍精度の参照場を float32 に丸めたまま測り、参照残差を汚染した事故がある
    # (`run_0027_s6_f64_from20`: `ro` 65193/65194 CV が変化、深部 max abs `roe` 1.168e-4)。
    # 指定しなかった場合はビット一致を求めず、**丸めで失われた大きさを報告**して続行する。
    bad, narrowed = [], []
    for n in moved:
        a = np.asarray(sv[n]); b = np.asarray(dv[n])
        if np.array_equal(a, b):
            continue
        if a.dtype.itemsize > b.dtype.itemsize:
            rel = np.max(np.abs(a.astype(np.float64) - b.astype(np.float64))) / max(np.max(np.abs(a)), 1e-300)
            narrowed.append((n, rel))
        else:
            bad.append(n)
    if bad:
        sys.exit(f"検査 NG: 写したのに SRC と一致しない: {bad}")
    if narrowed:
        print(f"⚠ 型の縮小 ({sv[moved[0]].dtype} -> {dv[moved[0]].dtype}) で丸めが入った。"
              f"**倍精度の場を種にするなら --keep-src-dtype を使うこと**:")
        for n, rel in narrowed:
            print(f"    {n:<10} 相対 {rel:.3e}")
    fsp.commit_inherit(d, species_plan)
    # 初期場の来歴: 座標の精度の引き上げなら記録し、旗なしの restart (同じ格子の場) なら前の記録を消す (docstring の 4)
    mpv.write_coord_upcast(d, upcast)
    print(f"species 属性  : {('継承 (species_input_unverified=%d)' % species_plan['species_input_unverified']) if species_plan else 'なし (未検証のまま)'}")
    print(f"VERDICT: OK ({len(moved)} 量を移した"
          f"{'、うち ' + str(len(narrowed)) + ' 量は型の縮小で丸めあり' if narrowed else '、SRC とビット一致'}"
          f"{'、' + mpv.coord_upcast_text(upcast) if upcast else ''})")
