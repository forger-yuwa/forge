#!/usr/bin/env python3
"""solverConfig.yaml の投入前チェック (procedures/recommended-settings.md の実体化; skill forge-config の手順 5)。

**残差を見ても気づけない設定ミス**を機械的に止めるのが目的。現在の検査:
  - `mesh.bndFirstOrder` は使用禁止 (AGENTS.md)
  - dual-time で `implicitRelax` < 1 (定常の推奨 0.7 を持ち込むと、残差は下がるのに遅いモードが収束せず
    同じ物理時刻の解が nSub 依存になる; §6 dual-time 節)
  - dual-time の `cfl_pseudo` が小さすぎる / `nSubIterDualTime` が少なすぎる (必要な nSub は cfl_pseudo で決まる)
  - 受動スカラの 2 次面移流 (`speciesFaceReconstruction` ≥ 2) が `convMethod 0` で不活性
  - FCT (`passiveFct 1`) が作動しない組み合わせ (SLAU 以外 / SFR < 2 / dual-time でない)
  - 凝縮 dual-time で `passiveScalarScheme 0` (モーメントに BDF 物理時間項が付かない)
  - `speciesFaceReconstruction` ≥ 2 と `speciesImplicitCoupling 0` の組 (定常で発散)
  - solver が読まないキー・**誤った節に書かれたキー** (WARN) と、起動時に拒否されるキー (FAIL)

使い方: check_solver_config.py RUN_DIR|solverConfig.yaml [...]   (VERDICT PASS/WARN/FAIL, exit 0/0/1)
WARN は「意図的ならよい」もの、FAIL は「まず直すべき」もの。
"""
import os
import sys


def load(path):
    import yaml
    p = os.path.join(path, 'solverConfig.yaml') if os.path.isdir(path) else path
    with open(p) as f:
        return p, (yaml.safe_load(f) or {})


_SPEC = None


def _spec():
    """solverConfig が受理する**完全修飾パス**の一覧を config_key_inventory から借りる。

    戻り (values, sections, rejected): `values` は値を読むパス、`sections` は節、`rejected` は
    起動時に落とすパス (削除済みキー・改名済みキー)。抽出に失敗したら None を返し、検査を飛ばす。"""
    global _SPEC
    if _SPEC is not None:
        return _SPEC
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location('cki', os.path.join(here, 'config_key_inventory.py'))
        cki = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cki)
        raw = open(cki.SRC_CPP, errors='replace').read() + '\n' + open(cki.SRC_HPP, errors='replace').read()
        src = cki.strip_comments(raw)
        keys = cki.extract_keys(src)
        gone = cki.removed_keys(src)
        allp = set(keys)
        sections = {p for p in allp if any(q != p and q.startswith(p + '.') for q in allp)}
        rejected = set(gone) | {p for p, d in keys.items() if d['kind'] == 'rejected'}
        values = (allp - sections - rejected)
        _SPEC = (values, sections, rejected)
    except Exception:
        _SPEC = (None, None, None)
    return _SPEC


def unknown_keys(y):
    """config に書かれた**完全修飾パス**を solverConfig の受理パスと突き合わせる (plan config-key-pruning §5.3 e)。

    戻り (fails, warns)。`fails` は solver が起動時に拒否するパス、`warns` は綴り違い・**誤配置** (正しい節が
    別にあるもの)・どこにも無いキー。末端名だけの照合では誤配置を拾えなかった (codex plan-2 M5)。"""
    values, sections, rejected = _spec()
    if values is None:
        return [], []
    by_leaf = {}
    for p in values:
        by_leaf.setdefault(p.split('.')[-1], []).append(p)
    fails, warns = [], []

    def walk(node, path):
        if not isinstance(node, dict):
            return
        for k, v in node.items():
            ks = str(k)
            p = '.'.join(path + [ks])
            if p in rejected:
                fails.append((p, 'solver が起動時に拒否するキー (削除済み・改名済み)。recommended-settings.md §9.1 の移行先を見ること'))
            elif p in values:
                continue                     # 値キー: 中身 (リスト・マップ) には立ち入らない
            elif p in sections:
                walk(v, path + [ks])
            elif ks in by_leaf:
                warns.append((p, '節の位置が違う (黙って無視される)。正しくは ' + ' / '.join(sorted(by_leaf[ks]))))
            else:
                warns.append((p, 'solverConfig が読まないキー (綴り違い・旧キー・別ブランチのキー)。黙って無視される'))

    walk(y, [])
    return fails, warns


def check(y):
    """戻り (fails, warns): それぞれ (キー, 説明) のリスト。"""
    fails, warns = [], []
    t = y.get('time', {}) or {}
    dT = t.get('deltaT', {}) or {}
    sp = y.get('space', {}) or {}
    cd = y.get('condensation', {}) or {}
    pp = y.get('physProp', {}) or {}
    mesh = y.get('mesh', {}) or {}
    solver = str(y.get('solver', ''))

    dual = (int(t.get('unsteady', 0) or 0) == 1 and int(t.get('dualTime', 0) or 0) == 1
            and int(t.get('timeIntegration', 0) or 0) == 11)
    nsub = int(t.get('nSubIterDualTime', 20) or 20)
    cflp = float(dT.get('cfl_pseudo', 0.0) or 0.0)
    relax = dT.get('implicitRelax')
    sfr = int(dT.get('speciesFaceReconstruction', 0) or 0)
    scheme = int(dT.get('passiveScalarScheme', 1) or 0)
    fct = int(dT.get('passiveFct', 1) or 0)
    conv = sp.get('convMethod')
    coupling = dT.get('speciesImplicitCoupling')

    ukf, ukw = unknown_keys(y)
    fails.extend(ukf)
    warns.extend(ukw)

    if mesh.get('bndFirstOrder') is not None:
        fails.append(('mesh.bndFirstOrder', '使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md)'))

    if dual:
        if relax is not None and float(relax) < 1.0:
            fails.append(('time.deltaT.implicitRelax',
                          f'dual-time で {relax} (定常の推奨値)。安定性は物理時間項が担うので効果が無く、遅いモードの収束だけ遅らせる。'
                          ' 残差ノルムには出ないまま解が nSub 依存になる (recommended-settings §6)。キーを消す (= 1.0) こと'))
        if cflp and cflp < 8.0:
            warns.append(('time.deltaT.cfl_pseudo',
                          f'{cflp:g} は小さく、必要な nSub が増える (cfl_pseudo 2 では nSub 40 でやっと収束、12 以上なら nSub 10 で足りる)。12–20 を推奨'))
        if nsub < 10:
            warns.append(('time.nSubIterDualTime', f'{nsub} は少ない。cfl_pseudo 12–20 と合わせて 10–20 を推奨'))
        if cd.get('condensation') and scheme != 1:
            fails.append(('time.deltaT.passiveScalarScheme',
                          '凝縮 dual-time で 0。モーメントに BDF 物理時間項が付かず、物理 Δt で積分される保証が無い (plan species-passive-scalar-unification §4.3)'))

    if sfr >= 2:
        if conv is not None and int(conv) == 0:
            fails.append(('space.convMethod',
                          f'speciesFaceReconstruction {sfr} なのに convMethod 0。面再構成が 1 次に落ちて S3 は不活性 (interp_dispatch)'))
        if coupling is not None and int(coupling) == 0:
            warns.append(('time.deltaT.speciesImplicitCoupling',
                          'S3 と coupling 0 の組は定常で組成せん断層から発散する (solver-settings.md)'))

    if scheme == 1 and fct == 1:
        why = []
        if solver not in ('SLAU', 'SLAU2'):
            why.append(f'solver {solver or "未指定"} (SLAU 系でない)')
        if sfr < 2:
            why.append(f'speciesFaceReconstruction {sfr} (< 2)')
        if not dual:
            why.append('dual-time でない')
        if why:
            warns.append(('time.deltaT.passiveFct', '1 だが作動しない: ' + ', '.join(why) + '。起動ログの `[passiveFct] active` で確認'))

    return fails, warns



def infos(y):
    """判定に影響しない注意 (INFO)。"""
    out = []
    mesh = y.get('mesh', {}) or {}
    sp = y.get('space', {}) or {}
    node = str(mesh.get('discretization', '')) == 'node'
    nwd = int(mesh.get('nodeWallDirichlet', 1) or 0) == 1
    slau = str(y.get('solver', '')).upper() in ('SLAU', 'SLAU2')
    if node and nwd and slau and 'slauWallNormalChi' not in sp:
        out.append(('space.slauWallNormalChi',
                    '省略 = auto で実効 1 (2026-09-26 に既定化)。2026-09-25 以前の結果を再現するには 0 を明記 '
                    '(plan convection-slau-wall-normal-chi-default、recommended-settings §1.0a)'))
    return out

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    if not args:
        print(__doc__)
        sys.exit(2)
    bad = False
    for a in args:
        try:
            p, y = load(a)
        except Exception as e:
            print(f'[{a}] cannot read solverConfig.yaml: {e}')
            bad = True
            continue
        fails, warns = check(y)
        print(f'--- {os.path.relpath(p)}')
        for k, m in fails:
            print(f'  FAIL {k}: {m}')
        for k, m in warns:
            print(f'  WARN {k}: {m}')
        for k, m in infos(y):
            print(f'  INFO {k}: {m}')
        if fails:
            bad = True
        print(f'  VERDICT: {"FAIL" if fails else ("WARN" if warns else "PASS")}')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
