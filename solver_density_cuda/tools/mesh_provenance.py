#!/usr/bin/env python3
"""変換済みメッシュの来歴 (変換時の実効条件) の記録と照合
(plan tooling-generator-defaults-and-launch-gate §4.3 の 2・§5.1 #7)。

`convertGmshToForge` は作業ディレクトリの `solverConfig.yaml`・`bcondConfig.yaml` を読んで HDF5 を書く。そのうち次は
**変換の結果を決め、格子に焼き込まれる**。`discretization`・`isAxisymmetric` 以外は forge の起動時の値が効かないので、run の config と
食い違っても気づかずに回る (`procedures/solver-settings.md` の「変換時」の節):

| 変換時の条件 | 効く先 (出典) |
| --- | --- |
| `mesh.discretization` | node なら双対を primal に置き換えて書き、primal を `/VIZMESH` に残す (`convertGmshToForge.cpp:44`) |
| `mesh.nodeInletCornerWall` | 入口の辺の壁節点側の半割面を壁の bcond へ (`gmshReader.hpp:1674`・`:2133`) |
| `mesh.axisCentroidShift` | 双対 CV の重心 (`gmshReader.hpp:2363`) |
| `mesh.renumber` | 節点の RCM 再番号 (`gmshReader.hpp:333`、`/MESH/RENUMBER_PERM`) |
| `mesh.isAxisymmetric` | 2D の node で `/PLANES/rSurfVect` を書く (`convertGmshToForge.cpp:50`) |
| `mesh.wallDistExtraPhysIDs` | `wall_dist` の壁点集合 (`calcWallDistance_kdtree.cpp:141`) |
| 境界の種類と physID | 壁 (`wall`・`wall_isothermal`) は `wall_dist` の壁点と入口角の帰属先 (`calcWallDistance_kdtree.cpp:142`、`gmshReader.hpp:1665`)、`inlet_*` は入口角の帰属元 (`gmshReader.hpp:1674`) |

**記録** `<h5>.provenance.json` (h5 の**実体** (realpath) のパスの末尾に `.provenance.json`。run が格子を symlink で置いても、
変換の時に実体の隣に書いた記録を引く。記録は実体の不変部分のハッシュに結び付くので、symlink の先や実体を差し替えれば照合で FAIL。
2026-10-15、plan §5.1 #23: 以前は symlink の名前の隣を探して、symlink の run を「記録の無いメッシュ」にしていた):

- `conversion`: 上の 6 キーの**要求値** (変換時の config に書いた値。書いていなければ `null` = 変換器の既定)。
  `mesh.isAxisymmetric` は solver と同じ優先順 (`mesh` > 旧 `physProp`) で読む。
- `bconds`: 変換時の bcondConfig の全境界 (physID・名前・種類・分類 `wall`/`inlet`/`other`)。
- `parts`: 変換成果物の**部分ごとのハッシュ** (`mesh_parts()`、版 `PARTS_SCHEME`)。`invariant` = 座標・接続・幾何・境界の所属・
  `wall_dist`・その他の合成。`/VALUE` (`wall_dist` 以外)・`/BCONDS/<id>/VALUE`・`/CHECKPOINT`・ルートの属性は**含めない**
  (restart・初期場の貼り付け・化学種の属性で変わる)。初期場は `initial_field` として別に記録する (照合には使わない)。
- `converter` (パス・sha256・終了コード)、`source` (`.msh` のパス・sha256)、`config_sha256` (変換時の 2 ファイル)、作成日時。

記録を書く前に、要求値と成果物の食い違い (`/VIZMESH` の有無と離散化、`/BCONDS/<id>` の `bcondKind` と bcondConfig、
`/MESH/RENUMBER_PERM` と `renumber`、軸対称でないのに `/PLANES/rSurfVect`) と、変換中の config の書き換えを調べ、あれば**書かない**
(変換器がその config を読んでいない)。変換器の版で決まる事実 (2D の node の軸対称で `/PLANES/rSurfVect` が有るか) は `artifact` に残す。
格子の種類の事実 (`node`・primal のセルの型ごとの数 `cell_types`、`mesh_kind()`) も `artifact` に残し、投入の関門が宣言の
プロファイルの適用条件 (次元・格子の種類・node/cell) と照合する (plan §5.1 #28。記録に `cell_types` が無ければ h5 から読む)。

**照合** (`check_run`、関門 `launch_gate.py` が呼ぶ予定。§5.1 #10) — FAIL:

1. run の `mesh.meshFileName` の h5 が無い・記録が無い・読めない・形式が違う (記録の無いメッシュ)。
2. 記録の `parts.invariant` が今の h5 と違う (別メッシュの記録・変換し直した h5 に旧い記録)。違う部分を名指しする。
3. 変換時の条件の要求値と run の config の要求値が違う。**省略は省略とだけ一致する** (変換器の既定と run の既定は
   バイナリで決まるので、今のソースの既定値で代用しない; plan §4.3 の 4 と同じ規則)。例外は `mesh.isAxisymmetric`:
   run が軸対称を要求 (1) するのに変換時に要求していなければ FAIL、逆 (軸対称で変換して平面で回す) は WARN、
   両方 1 でも 2D の node の格子に `/PLANES/rSurfVect` が無ければ (旧い変換器) WARN (forge は r̄_f·S_f に戻る)。
4. 境界: 変換時の physID が run の bcondConfig に無い、壁の分類 (`wall`・`wall_isothermal` かどうか) が変わった、
   入口の分類 (`inlet_*` かどうか) が変わった (変換時の `nodeInletCornerWall` が 0 と明記された node か cell のときは
   変換に効かないので WARN)。同じ分類の中の種類の違い (`inlet_Pressure` ↔ `inlet_uniformVelocity` など) は起動時の
   bcondConfig が効くので注記だけ。
5. `mesh.valueFileName` が別ファイルのとき: その `/VALUE/wall_dist` がメッシュの `wall_dist` と違う (長さ、または値。
   値の file の型が狭ければメッシュの値をその型に丸めて比べる; forge の res は読み込んだ型で書くため)。

使い方:
    mesh_provenance.py convert MSH OUT.h5 [--cwd DIR] [--converter BIN]   変換して <OUT.h5>.provenance.json を書く
    mesh_provenance.py check RUN_DIR [--config F] [--bcond F]              run と照合して VERDICT (PASS/WARN 0、FAIL 1)
    mesh_provenance.py parts H5                                             部分ごとのハッシュを出す
"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)

RECORD_FORMAT = 'forge-mesh-provenance'
RECORD_VERSION = 1
PARTS_SCHEME = 'mesh-parts-v1'
SUFFIX = '.provenance.json'

# 変換時の条件 (要求値として記録・照合するキー)。isAxisymmetric は別の規則 (上の 3)。
CONV_KEYS = ('mesh.discretization', 'mesh.nodeInletCornerWall', 'mesh.axisCentroidShift', 'mesh.renumber',
             'mesh.wallDistExtraPhysIDs', 'mesh.isAxisymmetric')
WALL_KINDS = ('wall', 'wall_isothermal')            # gmshReader.hpp:1665・calcWallDistance_kdtree.cpp:142

# 部分 (mesh_parts)。invariant はこの順の合成。
INVARIANT_PARTS = ('coord', 'topology', 'geometry', 'boundary', 'wall_dist', 'other')
_TOPOLOGY = {'MESH/CONNE', 'MESH/RENUMBER_PERM', 'VIZMESH/CONNE', 'PLANES/STRUCT', 'CELLS/STRUCT', 'CELLS/regionId',
             'DUAL/faceCells'}
_GEOMETRY = {'PLANES/surfVect', 'PLANES/surfArea', 'PLANES/centCoords', 'PLANES/rSurfVect', 'CELLS/volume',
             'CELLS/centCoords', 'DUAL/volume', 'DUAL/faceVect', 'DUAL/faceArea', 'DUAL/faceCent', 'DUAL/bnodeVect'}
_BOUNDARY_DUAL = {'DUAL/bnodeId', 'DUAL/bcondOffset', 'DUAL/bcondPhysID'}
_TOPOLOGY_GROUP_ATTRS = {'MESH', 'DUAL', 'VIZMESH'}
_CHUNK = 1 << 22                                    # 1 回に読む要素数の目安 (大きい格子でも全体を一度に持たない)


class ProvenanceError(Exception):
    """記録を作れない (変換の失敗・config と成果物の食い違い・変換中の書き換え・読めない入力)。"""


# ------------------------------------------------------------------------------------------------ 部分ごとのハッシュ

def part_of(path, kind, attr=None):
    """h5 の中の位置 (先頭の / なし) → 部分の名前。None = ハッシュに含めない。kind は 'dataset' か 'group'。
    attr は属性の名前 (属性のときだけ)。"""
    top = path.split('/', 1)[0] if path else ''
    if path == '':
        return None                                 # ルートの属性 (化学種の species_hash など。forge_species が書き換える)
    if top == 'VALUE':
        if path == 'VALUE/wall_dist' and attr is None and kind == 'dataset':
            return 'wall_dist'
        if path == 'VALUE' or path == 'VALUE/wall_dist':
            return None                             # /VALUE の属性・wall_dist の属性は持たない
        return 'initial_field'
    if top == 'CHECKPOINT':
        return 'initial_field'                      # dual-time の履歴 (状態)
    if top == 'BCONDS':
        segs = path.split('/')
        if len(segs) >= 3 and segs[2] == 'VALUE':
            return 'initial_field'                  # 境界の値 (変換時の bcondConfig の値。起動時は bcondConfig が効く)
        return 'boundary'                           # iPlanes・iBPlanes・iCells・vizBface*・属性 bcondKind
    if path == 'MESH/COORD' and attr is None:
        return 'coord'
    if path in _BOUNDARY_DUAL:
        return 'boundary'
    if path in _TOPOLOGY or (kind == 'group' and path in _TOPOLOGY_GROUP_ATTRS):
        return 'topology'
    if path in _GEOMETRY:
        return 'geometry'
    return 'other'                                  # 知らない位置も不変側に入れる (黙って外さない)


def _text(v):
    if isinstance(v, (bytes, np.bytes_)):
        try:
            return bytes(v).decode('utf-8')
        except UnicodeDecodeError:
            return None
    if isinstance(v, str):
        return v
    return None


def _canon_value(v):
    """属性・文字列のデータセットの値 → 決まったバイト列 (str と bytes は同じ文字列なら同じ)。"""
    t = _text(v)
    if t is not None:
        return b's:' + t.encode('utf-8')
    a = np.asarray(v)
    if a.dtype.kind in 'OSU':
        parts = []
        for x in a.ravel().tolist():
            tx = _text(x)
            parts.append(tx.encode('utf-8') if tx is not None else repr(x).encode())
        return b'S:' + repr(a.shape).encode() + b':' + b'\x00'.join(parts)
    a = np.ascontiguousarray(a)
    return b'n:' + a.dtype.str.encode() + repr(a.shape).encode() + b':' + a.tobytes()


def _update_dataset(h, ds):
    shape = ds.shape
    if shape is None:                               # null dataspace
        h.update(b'null')
        return
    if ds.dtype.kind in 'OSUV' or len(shape) == 0:
        h.update(_canon_value(ds[()]))
        return
    n = shape[0]
    row = int(np.prod(shape[1:])) if len(shape) > 1 else 1
    step = max(1, _CHUNK // max(1, row))
    for i in range(0, n, step):
        h.update(np.ascontiguousarray(ds[i:i + step]).tobytes())


def mesh_parts(h5):
    """変換成果物の部分ごとのハッシュ。h5 はパスか開いた h5py.File。戻り dict:

    - `scheme`: 版 (`PARTS_SCHEME`)
    - `parts`: {部分: sha256} — `coord` (/MESH/COORD を格納の型のまま)、`topology` (接続・CV の順・個数の属性)、
      `geometry` (面ベクトル・面積・重心・体積)、`boundary` (/BCONDS の所属と bcondKind、/DUAL の境界節点)、
      `wall_dist` (/VALUE/wall_dist の値)、`other` (知らない位置)、`initial_field` (照合に使わない初期場)
    - `invariant`: `INVARIANT_PARTS` の合成の sha256 (restart で変わらない部分)
    - `members`: {部分: [位置, ...]} (属性は `位置@名前`)
    - `counts`: /MESH・/DUAL の個数の属性

    #20 (`restart_field.py` の同一性の検査) が座標・接続・CV の順・境界の所属を部分ごとに比べられるよう、部分を分けて返す。"""
    import h5py
    own = not isinstance(h5, h5py.File)
    f = h5py.File(h5, 'r') if own else h5
    try:
        items = []                                  # (位置, 'dataset'/'group', obj)
        f.visititems(lambda name, obj: items.append(
            (name, 'dataset' if isinstance(obj, h5py.Dataset) else 'group', obj)))
        entries = {}                                # 部分 → [(鍵, 書き込み関数)]
        members = {}

        def add(part, key, fn):
            if part is None:
                return
            entries.setdefault(part, []).append((key, fn))
            members.setdefault(part, []).append(key)

        for name, kind, obj in items:
            if kind == 'dataset':
                add(part_of(name, kind), name,
                    (lambda o, nm: lambda h: (h.update(f'D|{nm}|{o.dtype.str}|{o.shape}|'.encode()),
                                              _update_dataset(h, o)))(obj, name))
            for an in sorted(obj.attrs.keys()):
                add(part_of(name, kind, an), f'{name}@{an}',
                    (lambda o, nm, a: lambda h: (h.update(f'A|{nm}@{a}|'.encode()),
                                                 h.update(_canon_value(o.attrs[a]))))(obj, name, an))
        parts = {}
        for part in INVARIANT_PARTS + ('initial_field',):
            h = hashlib.sha256()
            h.update(f'{PARTS_SCHEME}|{part}|'.encode())
            for key, fn in sorted(entries.get(part, []), key=lambda e: e[0]):
                fn(h)
            parts[part] = h.hexdigest()
        inv = hashlib.sha256(json.dumps({p: parts[p] for p in INVARIANT_PARTS}, sort_keys=True).encode()).hexdigest()
        counts = {}
        for g in ('MESH', 'DUAL'):
            if g in f:
                for an, v in f[g].attrs.items():
                    try:
                        counts[f'{g}.{an}'] = int(np.asarray(v))
                    except (TypeError, ValueError):
                        pass
        return {'scheme': PARTS_SCHEME, 'parts': parts, 'invariant': inv,
                'members': {k: sorted(v) for k, v in members.items()}, 'counts': counts}
    finally:
        if own:
            f.close()


# ------------------------------------------------------------------------------------------------ 初期場の来歴 (座標の精度の引き上げ)

# restart_field.py --allow-coord-upcast が DST の h5 のルートに書く記録 (JSON; plan §5.1 #20b)。ルートの属性は mesh_parts のハッシュに
# 入らない (part_of の '' → None) ので、記録を書いても格子の不変部分・来歴の照合は変わらない。投入の関門 (launch_gate.py) が
# LAUNCH_GATE.txt と RUN_PROVENANCE.txt に写す。
COORD_UPCAST_ATTR = 'restart_coord_upcast'


def coord_upcast_text(rec):
    """記録 → 1 行の説明「座標の精度の引き上げ (旧 → 新)」(旧・新 = 格子の不変部分の sha256 の先頭 16 桁)。"""
    return (f'座標の精度の引き上げ ({rec.get("src_coord_dtype")} の格子 {rec.get("src_invariant", "?")[:16]} → '
            f'{rec.get("dst_coord_dtype")} の格子 {rec.get("dst_invariant", "?")[:16]}; 格子の不変部分の sha256 '
            f'{rec.get("parts_scheme", PARTS_SCHEME)}; restart_field.py --allow-coord-upcast {rec.get("date", "")})')


def write_coord_upcast(h5, rec):
    """開いた h5py.File (書き込み可) のルートに記録を書く。rec が None なら前の記録を消す (同じ格子の場で上書きした初期場)。"""
    if rec is None:
        if COORD_UPCAST_ATTR in h5.attrs:
            del h5.attrs[COORD_UPCAST_ATTR]
        return
    h5.attrs[COORD_UPCAST_ATTR] = json.dumps(rec, ensure_ascii=False, sort_keys=True)


def read_coord_upcast(h5_path):
    """h5 のルートの記録 → dict (無ければ None)。読めない・形式が違えば ValueError。"""
    import h5py
    with h5py.File(h5_path, 'r') as f:
        if COORD_UPCAST_ATTR not in f.attrs:
            return None
        v = f.attrs[COORD_UPCAST_ATTR]
    t = _text(v)
    try:
        rec = json.loads(t) if t is not None else None
    except ValueError as e:
        raise ValueError(f'{COORD_UPCAST_ATTR} を JSON として読めない: {e}') from None
    if not isinstance(rec, dict) or rec.get('kind') != 'coord_upcast':
        raise ValueError(f'{COORD_UPCAST_ATTR} の形式が違う ({str(t)[:80]!r})')
    return rec


def provenance_path(h5_path):
    """記録のパス: メッシュの実体 (realpath) の末尾に `.provenance.json` (symlink の名前の隣ではない; 上の「記録」)。"""
    return os.path.realpath(str(h5_path)) + SUFFIX


def _mesh_desc(mesh_h5, name=None):
    """メッセージ用のメッシュの名前。symlink なら実体も書く。"""
    name = name or os.path.basename(mesh_h5)
    real = os.path.realpath(mesh_h5)
    return name if real == os.path.abspath(mesh_h5) else f'{name} (実体 {real} への symlink)'


# ------------------------------------------------------------------------------------------------ config・bcond の読み取り

def _load_yaml_strict(text, what):
    import yaml_strict
    try:
        doc = yaml_strict.load(text)
    except Exception as e:                          # 重複キー・merge key・構文の誤り
        raise ProvenanceError(f'{what} を読めない (検査できない入力): {e}') from None
    if not isinstance(doc, dict):
        raise ProvenanceError(f'{what} が YAML の写像でない (検査できない入力)')
    return doc


def _norm_int(v, key):
    if isinstance(v, bool):
        raise ProvenanceError(f'{key} が真偽値 ({v!r}); 整数で書く')
    if isinstance(v, int):
        return v
    if isinstance(v, float) and v.is_integer():
        return int(v)
    if isinstance(v, str):
        try:
            return int(v.strip().strip('"\''))
        except ValueError:
            pass
    raise ProvenanceError(f'{key} が整数でない ({v!r})')


def requested_conversion(cfg_doc):
    """solverConfig の写像 → {キー: 要求値 (書いていなければ None)}。書式の誤り (整数でない等) は ProvenanceError。"""
    mesh = cfg_doc.get('mesh')
    if mesh is not None and not isinstance(mesh, dict):
        raise ProvenanceError('solverConfig の mesh が写像でない')
    mesh = mesh or {}
    out = {}
    for key in CONV_KEYS:
        sub = key.split('.', 1)[1]
        if sub not in mesh:
            out[key] = None
            continue
        v = mesh[sub]
        if v is None:
            raise ProvenanceError(f'{key} が空 (`~`)。solver・変換器は空の値を読めない')
        if sub in ('discretization', 'renumber'):
            if not isinstance(v, str):
                raise ProvenanceError(f'{key} が文字列でない ({v!r})')
            out[key] = v.strip()                     # 変換器は大文字小文字を区別する (solverConfig.cpp:190・convertGmshToForge.cpp:33)
        elif sub == 'wallDistExtraPhysIDs':
            if not isinstance(v, (list, tuple)):
                raise ProvenanceError(f'{key} が列でない ({v!r})')
            out[key] = sorted({_norm_int(x, key) for x in v})       # 変換器は std::find で引くので順序・重複は効かない
        else:
            out[key] = _norm_int(v, key)
    if out['mesh.isAxisymmetric'] is None:          # solver と同じ: mesh 優先、無ければ旧 physProp (solverConfig.cpp:945・:953)
        pp = cfg_doc.get('physProp')
        if isinstance(pp, dict) and 'isAxisymmetric' in pp and pp['isAxisymmetric'] is not None:
            out['mesh.isAxisymmetric'] = _norm_int(pp['isAxisymmetric'], 'physProp.isAxisymmetric')
    return out


def bcond_class(kind):
    if kind in WALL_KINDS:
        return 'wall'
    if kind.startswith('inlet_'):                   # gmshReader.hpp:1674 (rfind("inlet_", 0) == 0)
        return 'inlet'
    return 'other'


def bcond_table(bc_doc):
    """bcondConfig の写像 → {physID: {'name', 'kind', 'class'}}。physID の重複・欠落は ProvenanceError。"""
    out = {}
    for name, d in bc_doc.items():
        if not isinstance(d, dict) or 'physID' not in d or 'kind' not in d:
            raise ProvenanceError(f'bcondConfig の {name!r} に physID か kind が無い (solver も読めない)')
        pid = _norm_int(d['physID'], f'{name}.physID')
        kind = d['kind']
        if not isinstance(kind, str) or not kind.strip():
            raise ProvenanceError(f'bcondConfig の {name!r} の kind が文字列でない ({kind!r})')
        if pid in out:
            raise ProvenanceError(f'bcondConfig の physID {pid} が 2 つの境界 ({out[pid]["name"]!r}・{name!r}) にある')
        kind = kind.strip()
        out[pid] = {'name': str(name), 'kind': kind, 'class': bcond_class(kind)}
    if not out:
        raise ProvenanceError('bcondConfig に境界が 1 つも無い')
    return out


def _sha_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ 成果物との突き合わせ

def _mesh_is_2d(f):
    """primal のセルが 2D の型か。変換器は gmsh の体積の有無で 2D/3D を決め、どちらかの次元のセルだけを書く
    (`gmshReader.hpp:324`・`:639`) ので、/VIZMESH/CONNE (node) か /MESH/CONNE (cell) の先頭のセルの型で見る
    (XDMF の 4 = 三角形・5 = 四角形、6〜9 = 3D)。読めなければ None。"""
    key = 'VIZMESH/CONNE' if 'VIZMESH/CONNE' in f else 'MESH/CONNE'
    if key not in f or f[key].shape is None or f[key].shape[0] == 0:
        return None
    code = int(f[key][0])
    return code in (4, 5) if code in (4, 5, 6, 7, 8, 9) else None


# primal のセルの型 (XDMF の Mixed の code → 名前・節点数。res_h5_to_vtu.XDMF2VTK と同じ)
XDMF_CELLS = {4: ('tri', 3), 5: ('quad', 4), 6: ('tet', 4), 7: ('pyramid', 5), 8: ('prism', 6), 9: ('hex', 8)}


def primal_cell_types(f):
    """primal のセルの型ごとの数 {名前: 数} (plan §5.1 #28: 適用条件の次元・格子の種類を成果物から取る)。node は /VIZMESH/CONNE
    (双対を primal に置き換えた後の primal; gmshReader.hpp:2659)、cell は /MESH/CONNE (XDMF の Mixed: code と節点の列)。
    知らない code (多角形・多面体)・列の長さの不整合・空は None (判定不能。一致扱いしない)。"""
    key = 'VIZMESH/CONNE' if 'VIZMESH/CONNE' in f else 'MESH/CONNE'
    if key not in f or f[key].shape is None:
        return None
    c = np.asarray(f[key])
    if c.ndim != 1 or c.size == 0 or int(c[0]) not in XDMF_CELLS:
        return None
    c0 = int(c[0])
    s = 1 + XDMF_CELLS[c0][1]
    # 1 種類だけの格子 (生産の格子はどれもこれ) は刻みで確かめる: 先頭から同じ型が続く限り次のセルは s の倍数の位置で始まるので、
    # 最初に違う型のセルは必ず s の倍数の位置の code で見つかる (取りこぼさない)
    if c.size % s == 0 and bool(np.all(c[::s] == c0)):
        return {XDMF_CELLS[c0][0]: int(c.size // s)}
    counts, i, n = {}, 0, int(c.size)
    while i < n:
        code = int(c[i])
        if code not in XDMF_CELLS:
            return None
        nm, nv = XDMF_CELLS[code]
        counts[nm] = counts.get(nm, 0) + 1
        i += 1 + nv
    return counts if i == n else None


def mesh_kind(h5_path):
    """変換の成果物 (h5) から格子の種類の事実: `node` (/VIZMESH の有無; artifact_problems と同じ印)・`cell_types`
    (primal_cell_types、判定不能は None)。来歴の記録の `artifact` にも同じ項目を書く (build_record)。"""
    import h5py
    with h5py.File(h5_path, 'r') as f:
        return {'node': 'VIZMESH' in f, 'cell_types': primal_cell_types(f)}


def artifact_problems(h5_path, conv, bconds):
    """要求値・境界の表と成果物の突き合わせ。戻り (食い違いの一覧, 成果物の事実 dict)。食い違いは「変換器がこの config・bcond を
    読んでいない」ことを示すもの (node の有無・境界の種類・再番号) で、あれば記録しない。変換器の版で決まるもの
    (`/PLANES/rSurfVect` の有無) は食い違いにせず事実として記録し、照合で WARN にする。"""
    import h5py
    probs = []
    with h5py.File(h5_path, 'r') as f:
        disc = conv['mesh.discretization'] or 'cell'      # 変換器は省略を cell と読む (solverConfig.cpp:194); 成果物で確かめる
        # node は双対を primal に置き換えて書き、primal のセルを /VIZMESH に退避する (gmshReader.hpp:2659)。/DUAL は置き換えの後
        # (dualBuilt = false, gmshReader.hpp:2487) には書かれないので、node の印は /VIZMESH で見る
        is_node = 'VIZMESH' in f
        if is_node != (disc == 'node'):
            probs.append(f'/VIZMESH が{"ある (node)" if is_node else "無い (cell)"}が、変換時の config の mesh.discretization は '
                         f'{conv["mesh.discretization"]!r} (変換器がこの config を読んでいない)')
        rn = conv['mesh.renumber']
        has_perm = 'MESH/RENUMBER_PERM' in f
        if has_perm != (rn == 'rcm'):
            probs.append(f'/MESH/RENUMBER_PERM が{"ある" if has_perm else "無い"}が、変換時の mesh.renumber は {rn!r}')
        has_rw = 'PLANES/rSurfVect' in f
        axi = conv['mesh.isAxisymmetric'] == 1
        if has_rw and not (axi and disc == 'node'):
            probs.append('/PLANES/rSurfVect があるが、変換時の config は node の軸対称でない')
        if 'BCONDS' not in f:
            probs.append('/BCONDS が無い')
        else:
            for k in f['BCONDS']:
                try:
                    pid = int(k)
                except ValueError:
                    probs.append(f'/BCONDS/{k} が physID の名前でない')
                    continue
                if pid not in bconds:
                    probs.append(f'/BCONDS/{pid} が変換時の bcondConfig に無い')
                    continue
                kind = _text(f['BCONDS'][k].attrs.get('bcondKind'))
                if kind != bconds[pid]['kind']:
                    probs.append(f'/BCONDS/{pid} の bcondKind {kind!r} が変換時の bcondConfig の {bconds[pid]["kind"]!r} と違う')
        if 'VALUE/wall_dist' not in f:
            probs.append('/VALUE/wall_dist が無い')
        facts = {'node': is_node, 'is2d': _mesh_is_2d(f), 'cell_types': primal_cell_types(f),
                 'rSurfVect': has_rw, 'renumber_perm': has_perm,
                 'coord_dtype': f['MESH/COORD'].dtype.str if 'MESH/COORD' in f else None,
                 'in_mesh': sorted(int(k) for k in f['BCONDS'] if str(k).isdigit()) if 'BCONDS' in f else []}
    return probs, facts


# ------------------------------------------------------------------------------------------------ 記録

def build_record(h5_path, cfg_text, bcond_text, converter=None, converter_rc=0, source=None, note=None):
    """変換の直後に呼ぶ: 変換時の config・bcond の文字列 (変換器を起動する前に読んだもの) と成果物から記録を作る。
    食い違いがあれば ProvenanceError (記録しない)。"""
    conv = requested_conversion(_load_yaml_strict(cfg_text, 'solverConfig.yaml'))
    bconds = bcond_table(_load_yaml_strict(bcond_text, 'bcondConfig.yaml'))
    probs, facts = artifact_problems(h5_path, conv, bconds)
    if probs:
        raise ProvenanceError('変換時の config と成果物が食い違う:\n    ' + '\n    '.join(probs))
    parts = mesh_parts(h5_path)
    in_mesh = set(facts.pop('in_mesh'))
    rec = {
        'format': RECORD_FORMAT, 'version': RECORD_VERSION,
        'mesh_file': os.path.basename(str(h5_path)),
        'created': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'tool': 'solver_density_cuda/tools/mesh_provenance.py',
        'conversion': conv,
        'bconds': [dict(physID=pid, in_mesh=(pid in in_mesh), **bconds[pid]) for pid in sorted(bconds)],
        'parts': {'scheme': parts['scheme'], 'invariant': parts['invariant'], **parts['parts']},
        'counts': parts['counts'],
        'artifact': facts,
        'config_sha256': {'solverConfig.yaml': hashlib.sha256(cfg_text.encode()).hexdigest(),
                          'bcondConfig.yaml': hashlib.sha256(bcond_text.encode()).hexdigest()},
        'converter': None, 'source': None,
    }
    if converter is not None:
        rec['converter'] = {'path': os.path.abspath(converter), 'returncode': int(converter_rc),
                            'sha256': _sha_file(converter) if os.path.isfile(converter) else None}
    if source is not None:
        rec['source'] = {'path': os.path.abspath(source), 'sha256': _sha_file(source) if os.path.isfile(source) else None}
    if note:
        rec['note'] = note
    return rec


def write_record(h5_path, rec):
    p = provenance_path(h5_path)
    tmp = p + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(rec, fh, indent=1, ensure_ascii=False)
        fh.write('\n')
    os.replace(tmp, p)
    return p


def default_converter():
    """convertGmshToForge の場所: FORGE_CONVERTER > FORGE_BIN と同じビルド > solver_density_cuda/build
    (design の `runner.converter_path()` と同じ順)。"""
    cand = os.environ.get('FORGE_CONVERTER')
    if cand:
        return cand
    fb = os.environ.get('FORGE_BIN')
    if fb and os.path.exists(os.path.join(os.path.dirname(fb), 'convertGmshToForge')):
        return os.path.join(os.path.dirname(fb), 'convertGmshToForge')
    return os.path.join(ROOT, 'solver_density_cuda', 'build', 'convertGmshToForge')


def _known_gpu_exit(r, out_h5):
    """AWS で既知の「HDF5 を書き終えた後の GPU の後片付けで落ちる」(design `runner.run_converter` と同じ判定)。"""
    if 'Write Input HDF5' not in (r.stdout or '') or 'GPUassert: invalid argument' not in (r.stderr or ''):
        return False
    try:
        import h5py
        with h5py.File(out_h5, 'r') as f:
            return 'MESH/COORD' in f and 'VALUE' in f
    except OSError:
        return False


def convert_and_record(msh, out_h5, cwd='.', converter=None, env=None):
    """変換器を cwd で回し (cwd の solverConfig.yaml・bcondConfig.yaml を読む)、記録を書く。戻り (記録のパス, 記録)。
    失敗は ProvenanceError (記録は書かず、旧い記録は消す)。"""
    cwd = os.path.abspath(cwd)
    converter = converter or default_converter()
    out_abs = out_h5 if os.path.isabs(out_h5) else os.path.join(cwd, out_h5)
    prov = provenance_path(out_abs)
    if os.path.exists(prov):
        os.remove(prov)                             # 成果物を作り直すので旧い記録は捨てる (失敗しても残さない)
    cfgp, bcp = os.path.join(cwd, 'solverConfig.yaml'), os.path.join(cwd, 'bcondConfig.yaml')
    for p in (cfgp, bcp):
        if not os.path.isfile(p):
            raise ProvenanceError(f'{p} が無い (変換器は作業ディレクトリのこのファイルを読む)')
    def read_bytes(p):
        with open(p, 'rb') as fh:
            return fh.read()
    cfg_bytes, bc_bytes = read_bytes(cfgp), read_bytes(bcp)
    try:
        cfg_text, bc_text = cfg_bytes.decode('utf-8'), bc_bytes.decode('utf-8')
    except UnicodeDecodeError as e:
        raise ProvenanceError(f'config が UTF-8 で読めない: {e}') from None
    # 変換の前に読めることを確かめる (読めない config で変換しても記録は作れない)
    requested_conversion(_load_yaml_strict(cfg_text, 'solverConfig.yaml'))
    bcond_table(_load_yaml_strict(bc_text, 'bcondConfig.yaml'))
    if not (os.path.isfile(converter) and os.access(converter, os.X_OK)):
        raise ProvenanceError(f'変換器 {converter} が無いか実行できない (--converter か FORGE_CONVERTER で指定)')
    # 変換器が書いたことを確かめる印 (旧い h5 が残っていて変換器が書かずに 0 で終わっても、旧い h5 を記録しない)
    before = os.stat(out_abs).st_mtime_ns if os.path.exists(out_abs) else None
    env = dict(os.environ if env is None else env)
    env['LD_LIBRARY_PATH'] = '/usr/lib/x86_64-linux-gnu/hdf5/serial:' + env.get('LD_LIBRARY_PATH', '')
    r = subprocess.run([converter, str(msh), str(out_h5)], cwd=cwd, env=env, capture_output=True, text=True)
    log = out_abs + '.convert.log'
    with open(log, 'w', encoding='utf-8') as fh:
        fh.write(f'# {converter} {msh} {out_h5} (cwd {cwd}) rc={r.returncode}\n--- stdout ---\n{r.stdout}\n--- stderr ---\n{r.stderr}')
    note = None
    if r.returncode != 0:
        if not _known_gpu_exit(r, out_abs):
            raise ProvenanceError(f'変換器が失敗 (rc={r.returncode}、ログ {log})')
        note = f'変換器は rc={r.returncode} で終わったが、終了時の既知の GPU エラーで HDF5 は書き終えている'
    if not os.path.isfile(out_abs):
        raise ProvenanceError(f'変換器は終わったが {out_abs} が無い (ログ {log})')
    if before is not None and os.stat(out_abs).st_mtime_ns == before:
        raise ProvenanceError(f'変換器は終わったが {out_abs} が書き換わっていない (旧い h5 のまま。ログ {log})')
    if read_bytes(cfgp) != cfg_bytes or read_bytes(bcp) != bc_bytes:
        raise ProvenanceError('変換の最中に solverConfig.yaml か bcondConfig.yaml が書き換わった (どちらで変換したか決まらない)')
    rec = build_record(out_abs, cfg_text, bc_text, converter=converter, converter_rc=r.returncode,
                       source=os.path.join(cwd, str(msh)) if not os.path.isabs(str(msh)) else str(msh), note=note)
    return write_record(out_abs, rec), rec


def load_record(path):
    with open(path, encoding='utf-8') as fh:
        rec = json.load(fh)
    if not isinstance(rec, dict) or rec.get('format') != RECORD_FORMAT:
        raise ProvenanceError(f'{path} は来歴の記録の形式でない (format {RECORD_FORMAT!r} が無い)')
    if rec.get('version') != RECORD_VERSION:
        raise ProvenanceError(f'{path} の版 {rec.get("version")!r} を読めない (この道具は {RECORD_VERSION})')
    parts = rec.get('parts')
    if not isinstance(parts, dict) or parts.get('scheme') != PARTS_SCHEME or not isinstance(parts.get('invariant'), str):
        raise ProvenanceError(f'{path} のハッシュの方式が {PARTS_SCHEME} でない (変換し直して記録を作る)')
    conv = rec.get('conversion')
    if not isinstance(conv, dict) or any(k not in conv for k in CONV_KEYS):
        raise ProvenanceError(f'{path} に変換時の条件 ({", ".join(CONV_KEYS)}) が揃っていない')
    bcs = rec.get('bconds')
    if not isinstance(bcs, list) or not bcs or not all(isinstance(b, dict) and {'physID', 'kind', 'class'} <= set(b) for b in bcs):
        raise ProvenanceError(f'{path} に境界の表が無い')
    return rec


# ------------------------------------------------------------------------------------------------ 照合

def _show(v):
    return '(書いていない)' if v is None else json.dumps(v, ensure_ascii=False)


def _same_file(a, b):
    try:
        return os.path.samefile(a, b)
    except OSError:
        return os.path.abspath(a) == os.path.abspath(b)


def check_value_wall_dist(mesh_h5, value_h5):
    """valueFileName が別ファイルのとき、その wall_dist がメッシュのものと同じか。戻り (問題の一覧, 注記の一覧)。
    値の file の型が狭いとき (float のビルドの res は float32) はメッシュの値をその型に丸めて比べる: forge は wall_dist を
    valueFileName から読む (variables.cpp:1023) ので、一致すればメッシュから直に読んだのと同じ値を float のビルドは使う。"""
    import h5py
    with h5py.File(mesh_h5, 'r') as fm, h5py.File(value_h5, 'r') as fv:
        if 'VALUE/wall_dist' not in fv:
            return [f'{os.path.basename(value_h5)} に /VALUE/wall_dist が無い (forge は wall_dist を valueFileName から読む: variables.cpp:1023)'], []
        a = np.asarray(fm['VALUE/wall_dist'])
        b = np.asarray(fv['VALUE/wall_dist'])
        if a.shape != b.shape:
            return [f'{os.path.basename(value_h5)} の wall_dist の長さ {b.shape} がメッシュの {a.shape} と違う (別のメッシュの場)'], []
        notes = []
        if a.dtype != b.dtype and a.dtype.kind == 'f' and b.dtype.kind == 'f' and b.dtype.itemsize < a.dtype.itemsize:
            ok = np.array_equal(a.astype(b.dtype), b)
            how = f'メッシュの {a.dtype} を {b.dtype} に丸めて比べた'
            notes.append(f'{os.path.basename(value_h5)} の wall_dist は {b.dtype} (メッシュは {a.dtype}): 倍精度のビルドでも丸めた壁距離で回る')
        else:
            ok = np.array_equal(a, b)
            how = '型のまま比べた'
        if not ok:
            nd = int(np.count_nonzero(a.astype(np.float64) != b.astype(np.float64)))
            return [f'{os.path.basename(value_h5)} の wall_dist がメッシュの wall_dist と {nd}/{a.size} 点で違う ({how})。'
                    '壁距離の条件 (壁の分類・wallDistExtraPhysIDs) が違う変換の場を使っている'], []
    return [], notes


def check_run(run_dir, cfg_path=None, bcond_path=None):
    """run の config・bcond とメッシュの来歴を照合する。戻り dict: fails/warns/notes [(キー, 説明)]、mesh・record・provenance のパス。"""
    fails, warns, notes = [], [], []
    out = {'fails': fails, 'warns': warns, 'notes': notes, 'mesh': None, 'provenance': None, 'record': None}
    cfg_path = cfg_path or os.path.join(run_dir, 'solverConfig.yaml')
    bcond_path = bcond_path or os.path.join(run_dir, 'bcondConfig.yaml')
    try:
        with open(cfg_path, encoding='utf-8') as fh:
            cfg = _load_yaml_strict(fh.read(), cfg_path)
        run_conv = requested_conversion(cfg)
    except (OSError, ProvenanceError) as e:
        fails.append(('(config)', f'run の config を検査できない: {e}'))
        return out
    mesh = cfg.get('mesh') or {}
    mname, vname = mesh.get('meshFileName'), mesh.get('valueFileName')
    if not isinstance(mname, str) or not mname.strip():
        fails.append(('mesh.meshFileName', '書いていない (検査できない)'))
        return out
    mesh_h5 = os.path.join(run_dir, mname.strip())
    out['mesh'] = mesh_h5
    if not os.path.isfile(mesh_h5):
        fails.append(('mesh.meshFileName', f'{mesh_h5} が無い'))
        return out
    prov = provenance_path(mesh_h5)
    out['provenance'] = prov
    if not os.path.isfile(prov):
        where = (os.path.basename(prov) if os.path.realpath(mesh_h5) == os.path.abspath(mesh_h5)
                 else f'{prov}; {_mesh_desc(mesh_h5, mname.strip())}、記録は実体の隣を探す')
        fails.append(('(provenance)', f'記録の無いメッシュ ({where} が無い)。'
                                      '`mesh_provenance.py convert` で変換し直す (今の config から過去の条件を推定して記録は作らない)'))
        return out
    try:
        rec = load_record(prov)
    except (OSError, ValueError, ProvenanceError) as e:
        fails.append(('(provenance)', f'記録を読めない: {e}'))
        return out
    out['record'] = rec
    try:
        cur = mesh_parts(mesh_h5)
    except Exception as e:
        fails.append(('(mesh)', f'{mesh_h5} を読めない (検査できない): {e}'))
        return out

    # --- 2. 成果物との結び付け ---
    if cur['invariant'] != rec['parts']['invariant']:
        diff = [p for p in INVARIANT_PARTS if rec['parts'].get(p) != cur['parts'][p]]
        hint = ''                                   # 手掛かり (判定は「違う」だけで決まる)
        if diff and set(diff) <= {'wall_dist'}:
            hint = ' — 違うのは wall_dist だけ: 壁距離の条件 (壁の分類・wallDistExtraPhysIDs) が違う変換の h5 に旧い記録か'
        elif 'coord' in diff:
            hint = ' — 座標が違う: 別のメッシュ (か renumber の違う変換) の記録か'
        elif diff:
            hint = (' — 座標は同じで接続・境界の所属などが違う: 同じ .msh を別の条件 (境界の分類・nodeInletCornerWall・'
                    'axisCentroidShift 等) で変換した h5 に旧い記録か、変換の後に書き換えた h5')
        fails.append(('(provenance)', f'記録が今の {_mesh_desc(mesh_h5)} と一致しない (違う部分: {", ".join(diff) or "合成のみ"})'
                                      f'{hint}。記録 {rec.get("mesh_file")} ({rec.get("created")})'))
    elif cur['parts']['initial_field'] != rec['parts'].get('initial_field'):
        notes.append(('(initial_field)', '初期場は変換の後に更新されている (restart・初期場の貼り付け。照合には使わない)'))

    # --- 3. 変換時の条件 ---
    conv = rec['conversion']
    for key in CONV_KEYS:
        c, r = conv.get(key), run_conv[key]
        if key == 'mesh.isAxisymmetric':
            c1, r1 = (c == 1), (r == 1)
            if r1 and not c1:
                fails.append((key, f'run は軸対称 (1) だが変換時は {_show(c)}: /PLANES/rSurfVect の無い格子では r̄_f·S_f のまま回る'
                                   ' (procedures/solver-settings.md の axisSegmentRWeight)。軸対称の config で変換し直す'))
            elif c1 and not r1:
                warns.append((key, f'軸対称 (1) で変換した格子を run は {_show(r)} で回す (変換で足した r 重みは使われない)'))
            elif c1 and r1:
                art = rec.get('artifact') or {}
                if art.get('node') and art.get('is2d') and not art.get('rSurfVect'):
                    warns.append((key, '2D の node の軸対称だが格子に /PLANES/rSurfVect が無い (旧い変換器): forge は r̄_f·S_f で回す'
                                       ' (procedures/solver-settings.md の axisSegmentRWeight)'))
            elif c != r:
                notes.append((key, f'変換時 {_show(c)}・run {_show(r)} (どちらも軸対称でない)'))
            continue
        if c != r:
            if c is None:
                why = '変換時は書いていない (変換器の既定。値を確かめられない)'
            elif r is None:
                why = f'変換時は {_show(c)} だが run は書いていない (run の既定はバイナリで決まる)'
            else:
                why = f'変換時は {_show(c)}、run は {_show(r)}'
            fails.append((key, f'{why}。変換時にだけ効く設定なので run の値は効かない — 同じ値を書くか、その config で変換し直す'))

    # --- 4. 境界 ---
    try:
        with open(bcond_path, encoding='utf-8') as fh:
            run_bc = bcond_table(_load_yaml_strict(fh.read(), bcond_path))
    except (OSError, ProvenanceError) as e:
        fails.append(('(bcond)', f'run の bcondConfig を検査できない: {e}'))
        run_bc = None
    if run_bc is not None:
        corner_off = (conv.get('mesh.discretization') != 'node') or (conv.get('mesh.nodeInletCornerWall') == 0)
        for b in rec['bconds']:
            pid = int(b['physID'])
            if not b.get('in_mesh', True):
                continue                            # 変換時の bcondConfig にだけあってメッシュに無い境界 (変換に効かない)
            tag = f'bcond.{pid}'
            if pid not in run_bc:
                fails.append((tag, f'変換時の境界 {b.get("name")!r} (physID {pid}, {b["kind"]}) が run の bcondConfig に無い'))
                continue
            rb = run_bc[pid]
            if rb['class'] == b['class']:
                if rb['kind'] != b['kind']:
                    notes.append((tag, f'種類 {b["kind"]} → {rb["kind"]} (同じ分類 {b["class"]}; 起動時は bcondConfig の種類が効く)'))
                continue
            if 'wall' in (rb['class'], b['class']):
                fails.append((tag, f'壁の分類が変わった ({b["kind"]} → {rb["kind"]}): wall_dist の壁点と入口角の帰属先は変換時の分類のまま'
                                   ' (calcWallDistance_kdtree.cpp:142・gmshReader.hpp:1665)。この bcond で変換し直す'))
            elif corner_off:
                warns.append((tag, f'入口の分類が変わった ({b["kind"]} → {rb["kind"]})。変換時は '
                                   f'{"cell" if conv.get("mesh.discretization") != "node" else "nodeInletCornerWall 0"} なので変換に効かない'
                                   ' (gmshReader.hpp:1674)'))
            else:
                fails.append((tag, f'入口の分類が変わった ({b["kind"]} → {rb["kind"]}): 変換時の nodeInletCornerWall '
                                   f'{_show(conv.get("mesh.nodeInletCornerWall"))} で入口角の半割面の帰属が決まっている'
                                   ' (gmshReader.hpp:1674)。この bcond で変換し直す'))
        extra = sorted(set(run_bc) - {int(b['physID']) for b in rec['bconds'] if b.get('in_mesh', True)})
        if extra:
            notes.append(('(bcond)', f'run の bcondConfig にだけある physID {extra} (メッシュに無いので使われない)'))

    # --- 5. valueFileName ---
    if not isinstance(vname, str) or not vname.strip():
        fails.append(('mesh.valueFileName', '書いていない (検査できない)'))
    else:
        value_h5 = os.path.join(run_dir, vname.strip())
        if not os.path.isfile(value_h5):
            fails.append(('mesh.valueFileName', f'{value_h5} が無い'))
        elif not _same_file(value_h5, mesh_h5):
            try:
                probs, vnotes = check_value_wall_dist(mesh_h5, value_h5)
                fails += [('mesh.valueFileName', p) for p in probs]
                notes += [('mesh.valueFileName', p) for p in vnotes]
            except Exception as e:
                fails.append(('mesh.valueFileName', f'{value_h5} を読めない (検査できない): {e}'))
    return out


def verdict(res):
    return 'FAIL' if res['fails'] else ('WARN' if res['warns'] else 'PASS')


def print_result(res, head):
    print(head)
    if res.get('mesh'):
        print(f'  メッシュ: {res["mesh"]}')
    rec = res.get('record')
    if rec:
        conv = rec['conversion']
        print(f'  記録    : {res["provenance"]} ({rec.get("created")}, 変換器 sha256 '
              f'{((rec.get("converter") or {}).get("sha256") or "不明")[:16]})')
        print('  変換時  : ' + ', '.join(f'{k.split(".", 1)[1]}={_show(conv[k])}' for k in CONV_KEYS))
    for k, m in res['fails']:
        print(f'  FAIL  {k}: {m}')
    for k, m in res['warns']:
        print(f'  WARN  {k}: {m}')
    for k, m in res['notes']:
        print(f'  注記  {k}: {m}')
    print(f'VERDICT: {verdict(res)}')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('convert', help='変換器を回して記録を書く')
    c.add_argument('msh')
    c.add_argument('out_h5')
    c.add_argument('--cwd', default='.', help='solverConfig.yaml・bcondConfig.yaml のある作業ディレクトリ (既定 .)')
    c.add_argument('--converter', help='convertGmshToForge (既定: FORGE_CONVERTER > FORGE_BIN のビルド > build/)')
    k = sub.add_parser('check', help='run の config・bcond と照合する')
    k.add_argument('run_dir')
    k.add_argument('--config')
    k.add_argument('--bcond')
    p = sub.add_parser('parts', help='部分ごとのハッシュを出す')
    p.add_argument('h5')
    a = ap.parse_args(argv)

    if a.cmd == 'convert':
        try:
            path, rec = convert_and_record(a.msh, a.out_h5, cwd=a.cwd, converter=a.converter)
        except ProvenanceError as e:
            print(f'[mesh_provenance] 記録しない: {e}')
            print('VERDICT: FAIL')
            return 1
        print(f'[mesh_provenance] 記録 {path}')
        print('  変換時: ' + ', '.join(f'{kk.split(".", 1)[1]}={_show(rec["conversion"][kk])}' for kk in CONV_KEYS))
        print(f'  境界  : ' + ', '.join(f'{b["physID"]}={b["kind"]}' for b in rec['bconds']))
        print(f'  不変部分 sha256: {rec["parts"]["invariant"]}')
        print('VERDICT: PASS')
        return 0
    if a.cmd == 'check':
        res = check_run(a.run_dir, a.config, a.bcond)
        print_result(res, f'[mesh_provenance] {a.run_dir}')
        return 1 if res['fails'] else 0
    try:
        r = mesh_parts(a.h5)
    except Exception as e:
        print(f'[mesh_provenance] {a.h5} を読めない: {e}')
        return 1
    print(json.dumps({'scheme': r['scheme'], 'invariant': r['invariant'], 'parts': r['parts'], 'counts': r['counts'],
                      'members': r['members']}, indent=1, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
