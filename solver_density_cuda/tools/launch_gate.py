#!/usr/bin/env python3
"""投入の関門: forge を起動する直前に、宣言・設定・メッシュの来歴・メッシュの品質・バイナリを照合する
(plan tooling-generator-defaults-and-launch-gate §4.3・§4.4・§4.7・§5.1 #10)。

`run_case.sh` が forge の起動の直前に `launch_gate.py <run_dir> --forge-bin <forge>` を呼び、終了コードが 0 でなければ
forge を起動しない。照合の中身はここに置く (`run_case.sh` は実行中に編集すると走行中の run を壊すので、呼び出しの数行だけにする)。

照合 (plan §4.3 の 1〜5):

1. **設定**: 宣言 (`<run_dir>/RECIPE.json`。`recipe_settings.py declare` か生成器が書く) があれば、宣言したプロファイル・段・逸脱で
   `recipe_settings.check_config` を回す (check_solver_config の検査、未知・誤配置のキー、YAML の重複キー・merge key、
   `fixed`・`explicit` の値と省略、`forbid`、プロファイルに無いキー、適用条件)。宣言の config の sha256 が今の
   `solverConfig.yaml` と違えば FAIL (宣言は別の config のもの。宣言し直す)。宣言が無い run は check_solver_config の検査だけを回す:
   その FAIL (削除・改名キー・`bndFirstOrder`・dual-time の緩和など、solver が起動を拒むかプロファイルの `forbid` と同じ類) と
   YAML の重複キー・merge key は FAIL、未知・誤配置のキーは宣言の欠落と同じ水準 (移行期 WARN・強制 FAIL)。
   **適用条件** (plan §5.1 #28): 宣言のプロファイルの `applies` のうち次元 (`2d_planar`・`2d_axisym`・`3d`)・格子の種類
   (`structured`・`unstructured`)・格子の node/cell は変換の成果物から取る (来歴の記録が今の h5 と一致すれば記録の `artifact` の
   `node`・`cell_types`、無ければ h5 の primal のセルの型; `gate_mesh_kind`)。導けない条件 (h5 が無い・読めない・型の混在など) は
   一致扱いせず、宣言した run では移行期でも FAIL (実験の宣言は NOTE)。
2. **メッシュの実効条件**: `mesh_provenance.check_run` (変換時の記録 `<h5>.provenance.json` と run の config・bcond の照合)。
   記録はメッシュの**実体** (realpath) の隣を探す (run が格子を symlink で置いても見つかる; plan §5.1 #23)。
3. **メッシュの品質**: `check_mesh_quality.py <h5> --bcond bcondConfig.yaml --ar-max 5000` (§4.4: AR 1000 を超えるセルは
   壁法線の構造層で条件を確かめられたときだけ 5000 まで、AR の外れ値に SOFT-PASS を当てない)。`bcondConfig.yaml` が無い・読めない
   (`--bcond` を渡せない) run は FAIL。結果の再利用は check_mesh_quality の鍵 (メッシュの不変部分・壁の physID と面・閾値・
   検査器の版) が全部一致したときだけ (`<h5>.quality.json`。実体の隣と run の中の名前の隣を読み、書くのは run の中だけ)。
4. **バイナリと省略値**: 起動するバイナリの sha256 と、宣言したプロファイルの内容のハッシュを記録する。このリポジトリの
   `solver_density_cuda/build/forge` でないバイナリ (対応の分からないバイナリ) は**拒否せず** WARN にする (旧・倍精度ビルドとの
   A/B を塞がない: plan §4.3 の 4)。省略の禁止は 1 の `fixed`・`explicit` の照合が担い、`key=~` の逸脱は実効値「未確定」で記録する。
5. 結果を `<run_dir>/LAUNCH_GATE.txt` (VERDICT・モード・宣言・プロファイル・段・逸脱・FAIL・WARN・メッシュのハッシュ・品質・
   バイナリ) に書き、`RUN_PROVENANCE.txt` に 1 行写す。初期場の h5 に `restart_field.py --allow-coord-upcast` の記録
   (座標の精度の引き上げ; plan §5.1 #20b) があれば両方に写す (判定は変えない)。

判定の規則 (検査できない入力はどのモードでも FAIL):

| 事象 | 移行期 (既定) | 強制 | `FORGE_LAUNCH_GATE=warn` |
| --- | --- | --- | --- |
| 宣言 (`RECIPE.json`) が無い | WARN | FAIL | 変えない |
| 宣言が読めない・形式が違う・理由の無い逸脱・宣言の config の sha256 が今の config と違う | FAIL | FAIL | 変えない |
| `fixed`・`explicit` の値の不一致・省略 (逸脱の宣言と config の値の食い違いを含む) | FAIL | FAIL | **WARN に落とす** |
| `forbid`・check_solver_config の FAIL・重複キー・merge key | FAIL | FAIL | 変えない |
| 未知/誤配置のキー・プロファイルに無いキー・適用条件の外 (宣言した run) | FAIL | FAIL | 変えない |
| 適用条件を照合できない (次元・格子の種類を成果物から取れない・知らない条件や値; 宣言した run) | FAIL | FAIL | 変えない |
| 未知/誤配置のキー (宣言の無い run) | WARN | FAIL | 変えない |
| `advisory` の範囲外・check_solver_config の WARN・宣言の後にプロファイルの内容が変わった | WARN | WARN | — |
| メッシュの記録 (`<h5>.provenance.json`) が無い | WARN | FAIL | 変えない |
| メッシュの来歴の照合の FAIL (記録の不一致・変換時の条件・境界の分類・valueFileName の wall_dist) | FAIL | FAIL | 変えない |
| 品質の FAIL / WARN・SOFT-PASS | FAIL / WARN | FAIL / WARN | 変えない |
| バイナリがこのリポジトリの build/forge でない・読めない | WARN (記録) | WARN (記録) | — |

**モード** (plan §4.7): 強制の開始の条件 (起動の経路の棚卸しの完了と、各経路の宣言の方法・段間の負例の試験が揃うこと) を
満たすまで、既定は「移行期」(宣言の無い run・記録の無いメッシュは WARN で起動する)。**切り替えは定数 `ENFORCE` の 1 か所**。
環境変数 `FORGE_LAUNCH_GATE` (`,` 区切り):

- `warn`: `fixed`・`explicit` の不一致・省略だけを WARN に落とす。`forbid`・宣言の欠落 (強制時)・メッシュの FAIL・検査不能は落とさない。
- `enforce`: その場で強制のモードにする (強制の開始の前に、経路ごとに先に確かめるため)。強制を弱める値は無い。

未知の値は FAIL (綴り違いで黙って既定に戻さない)。

使い方:
    launch_gate.py RUN_DIR [--forge-bin BIN] [--dry-run]
      --forge-bin  起動するバイナリ (既定: 環境変数 FORGE_BIN、無ければ solver_density_cuda/build/forge。run_case.sh と同じ)
      --dry-run    判定を出すだけで LAUNCH_GATE.txt・RUN_PROVENANCE.txt・品質の結果の再利用ファイルを書かない
終了コード: 0 = PASS・WARN (起動してよい)、1 = FAIL (起動しない)、2 = 引数の誤り。
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)

GATE_VERSION = 'launch-gate-v1'

# 強制の開始 (plan §4.7・§5.1 #3): 起動の経路の棚卸しが終わり、各経路の宣言の方法と段間の負例の試験が揃ったら True にする。
# **強制の既定を決めるのはここだけ**。False の間は宣言の無い run・記録の無いメッシュを WARN にして起動する。
ENFORCE = False

ENV = 'FORGE_LAUNCH_GATE'
MODE_TOKENS = ('warn', 'enforce')
RECIPE_FILE = 'RECIPE.json'
OUT_FILE = 'LAUNCH_GATE.txt'
PROV_FILE = 'RUN_PROVENANCE.txt'
CONFIG_FILE = 'solverConfig.yaml'
BCOND_FILE = 'bcondConfig.yaml'
DEFAULT_FORGE = os.path.join(ROOT, 'solver_density_cuda', 'build', 'forge')
QUALITY_TOOL = os.path.join(HERE, 'check_mesh_quality.py')

AREAS = ('モード', '宣言', '設定', '来歴', '品質', 'バイナリ', '関門')


def resolve_mode(environ=None):
    """環境変数からモードを決める。戻り ({'enforce', 'warn', 'env'}, 未知の値の一覧)。"""
    environ = os.environ if environ is None else environ
    raw = environ.get(ENV)
    toks = [t for t in re.split(r'[,\s]+', raw.strip().lower()) if t] if raw else []
    bad = [t for t in toks if t not in MODE_TOKENS]
    return {'enforce': bool(ENFORCE) or 'enforce' in toks, 'warn': 'warn' in toks, 'env': raw}, bad


def mode_text(mode):
    base = '強制' if mode['enforce'] else '移行期 (宣言の無い run・記録の無いメッシュは WARN)'
    if mode['enforce'] and not ENFORCE:
        base += ' (FORGE_LAUNCH_GATE=enforce)'
    if mode['warn']:
        base += '; FORGE_LAUNCH_GATE=warn (fixed・explicit の不一致・省略を WARN に落とす)'
    return base


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


class Findings:
    """判定の項目 (水準 FAIL/WARN/NOTE、区分、キー、説明)。"""

    def __init__(self):
        self.items = []

    def add(self, level, area, key, msg):
        assert level in ('FAIL', 'WARN', 'NOTE') and area in AREAS
        self.items.append((level, area, str(key), str(msg)))

    def count(self, level):
        return sum(1 for it in self.items if it[0] == level)

    def verdict(self):
        return 'FAIL' if self.count('FAIL') else ('WARN' if self.count('WARN') else 'PASS')


# ------------------------------------------------------------------------------------------------ 1. 設定と宣言

def _settings_undeclared(doc, bcond_doc, mode, out):
    """宣言の無い run (と読めない宣言) の設定の検査。check_solver_config の FAIL (削除・改名キー・bndFirstOrder・dual-time の緩和など、
    プロファイルの forbid と同じ類で solver が起動を拒むものを含む) は FAIL。未知・誤配置のキーは、宣言が無いあいだは
    「宣言の欠落」と同じ扱い (移行期 WARN・強制 FAIL)。宣言した run では recipe_settings.check_config が FAIL にする。"""
    import check_solver_config as csc
    if csc._spec()[0] is None:
        out.add('FAIL', '設定', '(受理キー)', 'solverConfig の受理キーの一覧を作れない (config_key_inventory の抽出に失敗; 検査できない)')
    _, uk_w = csc.unknown_keys(doc)
    cf, cw = csc.check(doc, bcond_doc)
    ukw = {k for k, _ in uk_w}
    for k, m in cf:
        out.add('FAIL', '設定', k, m)
    for k, m in uk_w:
        out.add('FAIL' if mode['enforce'] else 'WARN', '設定', k,
                m + (' (宣言した run と強制のモードでは FAIL)' if not mode['enforce'] else ''))
    for k, m in cw:
        if k not in ukw:
            out.add('WARN', '設定', k, m)


def _read_declaration(path):
    """RECIPE.json を読み、形を確かめる。戻り (宣言 dict, 問題の一覧)。"""
    import recipe_settings as RS
    try:
        with open(path, encoding='utf-8') as fh:
            decl = json.load(fh)
    except (OSError, ValueError) as e:
        return None, [f'{RECIPE_FILE} を読めない (検査できない): {e}']
    probs = []
    if not isinstance(decl, dict):
        return None, [f'{RECIPE_FILE} が JSON の写像でない']
    if decl.get('format') != RS.DECL_FORMAT or decl.get('version') != RS.DECL_VERSION:
        probs.append(f'形式 {decl.get("format")!r} 版 {decl.get("version")!r} を読めない '
                     f'(この関門は {RS.DECL_FORMAT} 版 {RS.DECL_VERSION})')
    kind = decl.get('kind')
    if kind not in ('profile', 'experiment'):
        probs.append(f'kind {kind!r} が profile / experiment でない')
    cfg = decl.get('config')
    if not isinstance(cfg, dict) or not isinstance(cfg.get('sha256'), str) or not cfg.get('sha256'):
        probs.append('config の sha256 が無い (どの config の宣言か分からない)')
    elif cfg.get('file') != CONFIG_FILE:
        probs.append(f'宣言の config は {cfg.get("file")!r} (forge が読むのは {CONFIG_FILE})')
    if kind == 'profile':
        for k in ('profile', 'stage'):
            if not isinstance(decl.get(k), str) or not decl.get(k).strip():
                probs.append(f'{k} が無い (プロファイルの宣言には profile と stage が要る)')
    if kind == 'experiment' and not str(decl.get('reason') or '').strip():
        probs.append('「実験」の宣言に理由 (reason) が無い')
    devs = decl.get('deviations', [])
    if not isinstance(devs, list):
        probs.append('deviations が列でない')
    else:
        for i, d in enumerate(devs):
            if not isinstance(d, dict) or not isinstance(d.get('key'), str) or not d.get('key'):
                probs.append(f'deviations[{i}] に key が無い')
            elif not str(d.get('reason') or '').strip():
                probs.append(f'逸脱 {d["key"]} に理由 (reason) が無い (理由の無い逸脱は宣言できない)')
    return decl, probs


def check_settings(run_dir, cfg_raw, doc, parse_fails, bcond_doc, mode, out, info, mesh_kind=None):
    """1. 設定と宣言。info に宣言・プロファイル・段・逸脱・プロファイルのハッシュを書く。mesh_kind: 格子の事実 (`gate_mesh_kind`)。
    宣言のプロファイルの適用条件のうち導けないもの (次元・格子の種類を成果物から取れない、など) は、宣言した run では移行期でも
    FAIL (判定不能を一致扱いしない; plan §5.1 #28)。`FORGE_LAUNCH_GATE=warn` でも落とさない (検査不能)。実験の宣言は NOTE。"""
    import recipe_settings as RS
    decl_path = os.path.join(run_dir, RECIPE_FILE)
    for k, m in parse_fails:
        out.add('FAIL', '設定', k, m)
    if not os.path.exists(decl_path):
        info['declaration'] = 'なし'
        out.add('FAIL' if mode['enforce'] else 'WARN', '宣言', RECIPE_FILE,
                '宣言が無い (推奨のプロファイルとの照合をしていない)。'
                '`python3 solver_density_cuda/tools/recipe_settings.py declare <run> --profile P --stage S` で宣言するか、'
                '推奨に従わない試し計算は `--experiment --reason "<1 行>"`'
                + ('' if mode['enforce'] else ' (強制の開始 (plan §4.7) までは WARN で起動する)'))
        if doc is not None:
            _settings_undeclared(doc, bcond_doc, mode, out)
        return
    decl, probs = _read_declaration(decl_path)
    info['declaration'] = f'{RECIPE_FILE} (sha256 {sha256_file(decl_path)})'
    for p in probs:
        out.add('FAIL', '宣言', RECIPE_FILE, p)
    if decl is None or probs:
        if doc is not None:
            _settings_undeclared(doc, bcond_doc, mode, out)
        return
    kind = decl['kind']
    info['declaration_kind'] = kind
    info['profile'] = decl.get('profile')
    info['stage'] = decl.get('stage')
    cur_sha = hashlib.sha256(cfg_raw).hexdigest() if cfg_raw is not None else None
    if cur_sha is not None and decl['config']['sha256'] != cur_sha:
        out.add('FAIL', '宣言', RECIPE_FILE,
                f'宣言の config の sha256 {decl["config"]["sha256"][:16]}… が今の {CONFIG_FILE} ({cur_sha[:16]}…) と違う: '
                '宣言の後に config が変わった (段の書き換えなど)。今の config で宣言し直す')
    if doc is None:
        return
    try:
        recipes = RS.load_recipes()
    except RS.RecipeError as e:
        out.add('FAIL', '設定', '(推奨設定)', f'推奨設定のファイルを読めない (検査できない): {e}')
        return
    info['recipes'] = f'{os.path.relpath(RS.RECIPES, ROOT)} (版 {recipes.get("version")}, sha256 {recipes.get("_sha256")})'
    devs = {d['key']: d.get('value') for d in decl.get('deviations') or []}
    info['deviations'] = list(decl.get('deviations') or [])
    import check_solver_config as csc
    if csc._spec()[0] is None:
        out.add('FAIL', '設定', '(受理キー)', 'solverConfig の受理キーの一覧を作れない (config_key_inventory の抽出に失敗; 検査できない)')
    try:
        r = RS.check_config(recipes, cfg_raw.decode('utf-8', errors='replace'), decl.get('profile'), decl.get('stage'), devs,
                            bcond_doc, experiment=(kind == 'experiment'), mesh=mesh_kind)
    except RS.RecipeError as e:
        out.add('FAIL', '設定', '(推奨設定)', f'プロファイルを解決できない (検査できない): {e}')
        return
    fixed = set(r.get('fixed_fails') or [])
    for k, m in r['fails']:
        if mode['warn'] and (k, m) in fixed:
            out.add('WARN', '設定', k, '[FORGE_LAUNCH_GATE=warn で FAIL から WARN に落とした] ' + m)
        else:
            out.add('FAIL', '設定', k, m)
    for k, m in r['warns']:
        out.add('WARN', '設定', k, m)
    for o in r.get('off_profile') or []:
        out.add('NOTE', '設定', o['key'], f'実験の宣言でプロファイルから外れた値 {RS.show(o["value"])} '
                                          f'(プロファイル {RS.show(o["profile_value"])})')
    for k, m in r.get('unchecked') or []:
        if kind == 'experiment':
            out.add('NOTE', '設定', k, f'実験の宣言: プロファイルの適用条件を照合していない ({m})')
        else:
            out.add('FAIL', '設定', k, f'プロファイルの適用条件を照合できない (判定不能は一致扱いしない; 移行期でも宣言した run は FAIL): {m}')
    res = r.get('profile')
    if res is not None:
        info['profile_sha256'] = res['sha256']
        declared = decl.get('profile_sha256')
        if not declared:
            out.add('WARN', '宣言', 'profile_sha256', '宣言にプロファイルの内容のハッシュが無い (宣言の時の内容と比べられない)')
        elif declared != res['sha256']:
            out.add('WARN', '宣言', 'profile_sha256',
                    f'宣言の後にプロファイルの内容が変わった (宣言 {declared[:16]}…、今 {res["sha256"][:16]}…)。照合は今の内容で行った')
    info['applied'] = r.get('applied') or []


# ------------------------------------------------------------------------------------------------ 2. メッシュの実効条件

def check_provenance(run_dir, mode, out, info):
    import mesh_provenance as mpv
    res = mpv.check_run(run_dir)
    info['mesh'] = res.get('mesh')
    rec = res.get('record')
    prov = res.get('provenance')
    no_record = rec is None and prov is not None and not os.path.isfile(prov) and res.get('mesh') and os.path.isfile(res['mesh'])
    for k, m in res['fails']:
        if no_record and k == '(provenance)':
            out.add('FAIL' if mode['enforce'] else 'WARN', '来歴', os.path.basename(prov),
                    m + ('' if mode['enforce'] else ' (強制の開始 (plan §4.7) までは WARN で起動する)'))
        else:
            out.add('FAIL', '来歴', k, m)
    for k, m in res['warns']:
        out.add('WARN', '来歴', k, m)
    for k, m in res['notes']:
        out.add('NOTE', '来歴', k, m)
    mesh = res.get('mesh')
    if mesh and os.path.isfile(mesh):
        mismatch = any(k == '(provenance)' for k, _ in res['fails'])
        if rec is not None and not mismatch:
            info['mesh_hash'] = f'{rec["parts"]["invariant"]} ({rec["parts"]["scheme"]}、来歴の記録と一致)'
        else:
            try:
                cur = mpv.mesh_parts(mesh)
                info['mesh_hash'] = (f'{cur["invariant"]} ({cur["scheme"]}、'
                                     + ('来歴の記録と不一致' if rec is not None else '来歴の記録なし') + ')')
            except Exception as e:                  # noqa: BLE001
                out.add('FAIL', '来歴', '(mesh)', f'{mesh} を読めない (検査できない): {e}')
    return res


def gate_mesh_kind(prov, info):
    """宣言のプロファイルの適用条件 (次元・格子の種類・格子の node/cell) と照合する格子の事実 (plan §5.1 #28)。来歴の記録が今の h5 と
    一致し、記録の成果物の事実 (`artifact`) に `node`・`cell_types` があればそれを使う。記録が無い・一致しない・記録に `cell_types`
    が無い (#28 より前の記録) ときは h5 そのものから読む (`mesh_provenance.mesh_kind`)。読めなければ None (判定不能: 宣言した run では
    照合できない適用条件として FAIL)。info['mesh_kind'] に事実と出所を書く。"""
    import mesh_provenance as mpv
    import recipe_settings as RS
    if not prov:
        return None
    mesh, rec = prov.get('mesh'), prov.get('record')
    consistent = rec is not None and not any(k == '(provenance)' for k, _ in prov.get('fails') or [])
    art = (rec.get('artifact') or {}) if consistent else {}
    if isinstance(art.get('node'), bool) and isinstance(art.get('cell_types'), dict):
        facts = {'node': art['node'], 'cell_types': art['cell_types']}
        info['mesh_kind'] = f'{RS.mesh_kind_text(facts)} (来歴の記録の成果物の事実)'
        return facts
    if not mesh or not os.path.isfile(mesh):
        info['mesh_kind'] = '不明 (メッシュの h5 が無い)'
        return None
    try:
        facts = mpv.mesh_kind(mesh)
    except Exception as e:                          # noqa: BLE001
        info['mesh_kind'] = f'不明 (h5 を読めない: {e})'
        return None
    why = ('来歴の記録に cell_types が無い' if consistent else
           ('来歴の記録と一致しない' if rec is not None else '来歴の記録が無い'))
    info['mesh_kind'] = f'{RS.mesh_kind_text(facts)} (h5 から読んだ; {why})'
    return facts


def note_initial_field(run_dir, doc, out, info):
    """初期場の来歴: run が読む初期場の h5 (`mesh.valueFileName`、無ければ `mesh.meshFileName`) に `restart_field.py
    --allow-coord-upcast` の記録 (`mesh_provenance.COORD_UPCAST_ATTR`) があれば LAUNCH_GATE.txt と RUN_PROVENANCE.txt に写す
    (plan §5.1 #20b。`run_case.sh` は起動ごとに RUN_PROVENANCE.txt を書き直すので、記録は h5 に置いてここで写す)。判定は変えない。"""
    import mesh_provenance as mpv
    mesh = (doc or {}).get('mesh') if isinstance(doc, dict) else None
    if not isinstance(mesh, dict):
        return
    name = mesh.get('valueFileName') or mesh.get('meshFileName')
    if not isinstance(name, str) or not name.strip():
        return
    path = os.path.join(run_dir, name.strip())
    if not os.path.isfile(path):
        return                                      # 無いファイルは来歴の照合 (check_provenance) が扱う
    try:
        rec = mpv.read_coord_upcast(path)
    except (OSError, ValueError) as e:
        out.add('NOTE', '初期場', os.path.basename(path), f'初期場の来歴の記録を読めない: {e}')
        return
    if rec is not None:
        info['initial_field'] = f'{os.path.basename(path)}: {mpv.coord_upcast_text(rec)}'


# ------------------------------------------------------------------------------------------------ 3. メッシュの品質

def check_quality(mesh, bcond_path, out, info, dry_run=False):
    """check_mesh_quality.py --bcond を子プロセスで回す (再利用は検査器の鍵に任せる)。"""
    import check_mesh_quality as cmq
    ar_max = cmq.AR_RELAX_CAP
    fd, js = tempfile.mkstemp(prefix='launch_gate_quality_', suffix='.json')
    os.close(fd)
    os.unlink(js)
    cmd = [sys.executable, QUALITY_TOOL, mesh, '--bcond', bcond_path, '--ar-max', '%g' % ar_max, '--json', js]
    if dry_run:
        cmd.append('--no-cache')
    info['quality_cmd'] = 'check_mesh_quality.py %s --bcond %s --ar-max %g%s' % (
        os.path.basename(mesh), os.path.basename(bcond_path), ar_max, ' --no-cache' if dry_run else '')
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
        lines = (p.stdout + p.stderr).rstrip().splitlines()
        info['quality_lines'] = lines
        summ = None
        if os.path.isfile(js):
            with open(js, encoding='utf-8') as fh:
                summ = json.load(fh)
    except Exception as e:                          # noqa: BLE001
        out.add('FAIL', '品質', '(mesh)', f'品質の検査を回せない (検査できない): {e}')
        return
    finally:
        if os.path.exists(js):
            os.unlink(js)
    tail = next((ln for ln in reversed(lines) if ln.startswith('VERDICT')), lines[-1] if lines else '(出力なし)')
    if not isinstance(summ, dict) or 'verdict' not in summ:
        out.add('FAIL', '品質', os.path.basename(mesh), f'品質の判定を得られない (検査できない; 終了 {p.returncode}): {tail}')
        return
    v = str(summ.get('verdict'))
    info['quality'] = f'{v} (再利用 {"あり" if summ.get("reused") else "なし"}; 終了 {p.returncode})'
    if (p.returncode == 0) != (v != 'FAIL'):
        out.add('FAIL', '品質', os.path.basename(mesh), f'品質の判定 {v} と終了コード {p.returncode} が食い違う (検査できない)')
    elif v == 'FAIL':
        out.add('FAIL', '品質', os.path.basename(mesh), tail)
    elif v in ('WARN', 'SOFT-PASS'):
        out.add('WARN', '品質', os.path.basename(mesh), tail)
    elif v != 'PASS':
        out.add('FAIL', '品質', os.path.basename(mesh), f'品質の判定 {v!r} を知らない (検査できない)')


# ------------------------------------------------------------------------------------------------ 4. バイナリ

def check_binary(forge_bin, out, info, environ=None):
    environ = os.environ if environ is None else environ
    path = forge_bin or environ.get('FORGE_BIN') or DEFAULT_FORGE
    info['forge_bin'] = path
    if not os.path.isfile(path):
        info['forge_sha256'] = '不明 (無い・読めない)'
        out.add('WARN', 'バイナリ', path, 'バイナリが無い・読めない (sha256 を記録できない)')
        return
    try:
        info['forge_sha256'] = sha256_file(path)
    except OSError as e:
        info['forge_sha256'] = '不明 (読めない)'
        out.add('WARN', 'バイナリ', path, f'バイナリを読めない (sha256 を記録できない): {e}')
        return
    try:
        info['forge_mtime'] = datetime.datetime.fromtimestamp(os.path.getmtime(path)).astimezone().isoformat(timespec='seconds')
    except OSError:
        pass
    if os.path.realpath(path) != os.path.realpath(DEFAULT_FORGE):
        out.add('WARN', 'バイナリ', path,
                f'対応の分からないバイナリ (このリポジトリの {os.path.relpath(DEFAULT_FORGE, ROOT)} でない)。省略したキーの値・既定は'
                'そのバイナリで決まる。拒否はしない (旧・倍精度ビルドとの A/B のため; plan §4.3 の 4)。sha256 を記録した')


# ------------------------------------------------------------------------------------------------ 本体

def run_gate(run_dir, forge_bin=None, environ=None, dry_run=False):
    """関門を回す。戻り (終了コード, Findings, info)。LAUNCH_GATE.txt と RUN_PROVENANCE.txt の 1 行を書く (dry_run 以外)。"""
    run_dir = os.path.abspath(run_dir)
    out = Findings()
    info = {'run_dir': run_dir, 'date': datetime.datetime.now().astimezone().isoformat(timespec='seconds')}
    mode, bad = resolve_mode(environ)
    info['mode'] = mode_text(mode)
    info['mode_flags'] = mode
    try:
        for t in bad:
            out.add('FAIL', 'モード', ENV, f'未知の値 {t!r} (使えるのは {", ".join(MODE_TOKENS)}; 綴り違いで黙って既定に戻さない)')
        if not os.path.isdir(run_dir):
            out.add('FAIL', '関門', run_dir, 'run ディレクトリが無い')
            return _finish(run_dir, out, info, dry_run)
        cfg_path = os.path.join(run_dir, CONFIG_FILE)
        cfg_raw, doc, parse_fails = None, None, []
        try:
            with open(cfg_path, 'rb') as fh:
                cfg_raw = fh.read()
        except OSError as e:
            parse_fails = [('(config)', f'{CONFIG_FILE} を読めない (検査できない): {e}')]
        if cfg_raw is not None:
            import recipe_settings as RS
            info['config_sha256'] = hashlib.sha256(cfg_raw).hexdigest()
            doc, parse_fails = RS.parse_config(cfg_raw.decode('utf-8', errors='replace'))
        bcond_path = os.path.join(run_dir, BCOND_FILE)
        bcond_doc = None
        if not os.path.isfile(bcond_path):
            out.add('FAIL', '品質', BCOND_FILE, f'{BCOND_FILE} が無い: 壁の bcond が分からず、品質の検査に --bcond を渡せない (検査できない)')
            bcond_path = None
        else:
            import yaml_strict
            try:
                with open(bcond_path, encoding='utf-8') as fh:
                    bcond_doc = yaml_strict.load(fh.read())
                if not isinstance(bcond_doc, dict):
                    raise ValueError('YAML の写像でない')
            except Exception as e:                  # noqa: BLE001
                out.add('FAIL', '品質', BCOND_FILE, f'{BCOND_FILE} を読めない (検査できない): {e}')
                bcond_doc, bcond_path = None, None
        # 2. メッシュの実効条件 (設定の照合より先に回し、宣言のプロファイルの適用条件 (次元・格子の種類) と照合する格子の事実を
        # 取る: plan §5.1 #28。判定の項目は従来どおり設定の後に並べる)
        pout = Findings()
        prov = check_provenance(run_dir, mode, pout, info) if doc is not None else None
        mesh_facts = gate_mesh_kind(prov, info)
        # 1. 設定と宣言
        check_settings(run_dir, cfg_raw, doc, parse_fails, bcond_doc, mode, out, info, mesh_kind=mesh_facts)
        out.items += pout.items
        note_initial_field(run_dir, doc, out, info)
        # 3. メッシュの品質
        mesh = (prov or {}).get('mesh')
        if mesh and os.path.isfile(mesh) and bcond_path is not None:
            check_quality(mesh, bcond_path, out, info, dry_run=dry_run)
        # 4. バイナリ
        check_binary(forge_bin, out, info, environ)
    except Exception as e:                          # noqa: BLE001
        out.add('FAIL', '関門', '(例外)', f'関門の内部の例外 (検査できない): {type(e).__name__}: {e}')
    return _finish(run_dir, out, info, dry_run)


def render(out, info):
    """LAUNCH_GATE.txt の本文。"""
    L = [f'# 投入の関門 ({GATE_VERSION}; solver_density_cuda/tools/launch_gate.py、plan tooling-generator-defaults-and-launch-gate §4.3)',
         f'date        : {info.get("date")}',
         f'run_dir     : {info.get("run_dir")}',
         f'mode        : {info.get("mode")}',
         f'declaration : {info.get("declaration", "-")}' + (f' (kind {info["declaration_kind"]})' if info.get('declaration_kind') else ''),
         f'profile     : {info.get("profile") or "-"}' + (f' (sha256 {info["profile_sha256"]})' if info.get('profile_sha256') else ''),
         f'stage       : {info.get("stage") or "-"}',
         f'config      : {CONFIG_FILE} (sha256 {info.get("config_sha256", "不明")})']
    if info.get('recipes'):
        L.append(f'recipes     : {info["recipes"]}')
    for d in info.get('deviations') or []:
        eff = ''
        for a in info.get('applied') or []:
            if a.get('key') == d.get('key') and a.get('effective'):
                eff = f'; 実効値 {a["effective"]}'
        val = '~ (書かない)' if d.get('value') is None else json.dumps(d.get('value'), ensure_ascii=False)
        L.append(f'deviation   : {d.get("key")} = {val} (プロファイル {json.dumps(d.get("profile_value"), ensure_ascii=False)}'
                 f'{eff}; 理由 {d.get("reason")}' + (f'; 根拠 {d["basis"]}' if d.get('basis') else '') + ')')
    L.append(f'mesh        : {info.get("mesh") or "-"}')
    L.append(f'mesh_hash   : {info.get("mesh_hash", "-")}')
    if info.get('mesh_kind'):
        L.append(f'mesh_kind   : {info["mesh_kind"]}')
    if info.get('initial_field'):
        L.append(f'initial_field: {info["initial_field"]}')
    L.append(f'quality     : {info.get("quality", "-")}' + (f' [{info["quality_cmd"]}]' if info.get('quality_cmd') else ''))
    L.append(f'forge_bin   : {info.get("forge_bin", "-")}')
    L.append(f'forge_sha256: {info.get("forge_sha256", "-")}')
    if info.get('forge_mtime'):
        L.append(f'forge_mtime : {info["forge_mtime"]}')
    for lv in ('FAIL', 'WARN', 'NOTE'):
        for level, area, key, msg in out.items:
            if level == lv:
                L.append(f'{lv:<4}        : [{area}] {key}: {msg}')
    q = info.get('quality_lines') or []
    if q:
        L.append('quality_out :')
        L += ['  | ' + ln for ln in q[-40:]]
    L.append(f'VERDICT: {out.verdict()} (FAIL {out.count("FAIL")}・WARN {out.count("WARN")}; '
             + ('forge を起動しない' if out.verdict() == 'FAIL' else '起動してよい') + ')')
    return '\n'.join(L) + '\n'


def _finish(run_dir, out, info, dry_run):
    text = render(out, info)
    sys.stdout.write(text)
    v = out.verdict()
    if not dry_run and os.path.isdir(run_dir):
        try:
            p = os.path.join(run_dir, OUT_FILE)
            tmp = f'{p}.tmp{os.getpid()}'
            with open(tmp, 'w', encoding='utf-8') as fh:
                fh.write(text)
            os.replace(tmp, p)
            mf = info.get('mode_flags') or {}
            line = (f'launch_gate : VERDICT {v} ({"強制" if mf.get("enforce") else "移行期"}'
                    f'{"; FORGE_LAUNCH_GATE=warn" if mf.get("warn") else ""}; '
                    f'profile {info.get("profile") or "-"}; stage {info.get("stage") or "-"}; '
                    f'FAIL {out.count("FAIL")}・WARN {out.count("WARN")}; forge_sha256 {str(info.get("forge_sha256", "-"))[:16]}; '
                    + (f'初期場 {info["initial_field"]}; ' if info.get('initial_field') else '')
                    + f'詳細 {OUT_FILE})\n')
            with open(os.path.join(run_dir, PROV_FILE), 'a', encoding='utf-8') as fh:
                fh.write(line)
        except OSError as e:
            print(f'[launch_gate] {OUT_FILE} / {PROV_FILE} を書けない: {e} -> FAIL')
            return 1, out, info
    return (1 if v == 'FAIL' else 0), out, info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('run_dir')
    ap.add_argument('--forge-bin', help='起動するバイナリ (既定: FORGE_BIN、無ければ solver_density_cuda/build/forge)')
    ap.add_argument('--dry-run', action='store_true', help='判定を出すだけでファイルを書かない')
    a = ap.parse_args(argv)
    rc, _, _ = run_gate(a.run_dir, a.forge_bin, dry_run=a.dry_run)
    return rc


if __name__ == '__main__':
    sys.exit(main())
