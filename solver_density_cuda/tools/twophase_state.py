#!/usr/bin/env python3
"""二相拡散 (`condensation.condTwoPhaseDiffusion`) の**実効状態**を solverConfig.yaml (+ bcondConfig.yaml) から判定する。

plans/active/condensation-two-phase-default.md §4-2 の判定表の Python 側。ソルバ側の正本は
`cuda_forge/condensationTransport_d.cuh` の `condTwoPhaseDiffusionClassify` (起動行 `[twophase]`・`res_*.h5` 属性
`twophase_diffusion_state` / `twophase_diffusion_effective`)。**判定を変えたら両方を変える** (試験
`tests/unit/test_twophase_state.py` が既定値の一致を見る)。

実効状態 (指定値は見ない。構成だけで決まる; 判定順は (a) → (b) → (c)):
  active        包絡内。指定 ON で作動する。
  inactive-a    物理が同一 (凝縮 OFF・viscMethod 0)。どの指定でもビット一致。dual-time + Euler もここ。
  inactive-b    モデルが構造的に適用できない (CPG carrier・pure 凝縮・condEquilibrium 2)。指定 ON なら毎回 WARNING。
  unsupported-c 実装が未対応 (dual-time・RK/陽解法・speciesImplicitCoupling 2・passiveScalarScheme 0・condEquilibrium 1・
                condLimiterMode 0・nCondSpecies ≥ 2・軸対称・周期)。指定 ON ならソルバはエラー終了。

使い方: twophase_state.py RUN_DIR|solverConfig.yaml   (判定結果を 1 行出す)
"""
import os
import sys

# 省略時の既定。ソルバの `solverConfig::kCondTwoPhaseDiffusionDefault` (input/solverConfig.hpp) と同じ値にする。
# S1-c (既定 ON) で両方を 1 にする。省略を ON と数えるか・明示 0 を WARN にするかはこの値で決まる。
DEFAULT = 0

STATES = ('active', 'inactive-a', 'inactive-b', 'unsupported-c')


def _int(d, k, default):
    v = (d or {}).get(k, default)
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def _has_periodic(bcond):
    """bcondConfig (dict) に kind: periodic の境界があるか。"""
    if not isinstance(bcond, dict):
        return False
    for v in bcond.values():
        if isinstance(v, dict) and str(v.get('kind', '')).strip() == 'periodic':
            return True
    return False


def classify(y, bcond=None, default=None):
    """戻り dict: requested ('omitted'/'0'/'1'), request_on (bool), state, effective (0/1), reason。

    y は solverConfig.yaml を読んだ dict、bcond は bcondConfig.yaml の dict (周期の判定に使う; 無ければ非周期とみなす)。
    default は省略時の既定 (None ならモジュールの DEFAULT)。"""
    if default is None:
        default = DEFAULT
    y = y or {}
    cd = y.get('condensation') or {}
    pp = y.get('physProp') or {}
    t = y.get('time') or {}
    dT = t.get('deltaT') or {}
    mesh = y.get('mesh') or {}

    given = isinstance(cd, dict) and 'condTwoPhaseDiffusion' in cd
    req_val = _int(cd, 'condTwoPhaseDiffusion', default) if given else default
    request_on = (req_val == 1)
    requested = str(req_val) if given else 'omitted'

    condensation = _int(cd, 'condensation', 0)
    n_cond = _int(cd, 'nCondSpecies', 0) if condensation == 1 else 0
    thermal = _int(pp, 'thermalMethod', 0)
    visc = _int(pp, 'viscMethod', 0)
    species = pp.get('species') or []
    names = []
    for s in species:
        names.append(str(s.get('name')) if isinstance(s, dict) else str(s))
    n_species = len(names)
    if thermal == 2 and n_species == 0:
        n_species = 1                               # ソルバの既定 (単成分 N2)
    gas = _int(cd, 'condGasSpecies', -1)
    if cd.get('condensationSpecies') is not None:   # 名前が正本 (大文字小文字無視)
        up = str(cd['condensationSpecies']).upper()
        gas = next((i for i, n in enumerate(names) if n.upper() == up), -1)
    try:
        vapor_mf = float(cd.get('condVaporMassFraction', -1.0))
    except (TypeError, ValueError):
        vapor_mf = -1.0
    equil = _int(cd, 'condEquilibrium', 0)
    lim_mode = _int(cd, 'condLimiterMode', 1)
    unsteady = _int(t, 'unsteady', 0)
    dual = _int(t, 'dualTime', 0)
    ti = _int(t, 'timeIntegration', 0)
    coupling = _int(dT, 'speciesImplicitCoupling', 0)
    scheme = _int(dT, 'passiveScalarScheme', 1)
    axisym = _int(mesh, 'isAxisymmetric', _int(pp, 'isAxisymmetric', 0))
    periodic = _has_periodic(bcond)

    def out(state, reason):
        return {'requested': requested, 'request_on': request_on, 'state': state,
                'effective': 1 if (request_on and state == 'active') else 0, 'reason': reason}

    # (a) 物理が同一
    if condensation != 1 or n_cond < 1:
        return out('inactive-a', 'condensation is off')
    if visc == 0:
        return out('inactive-a', 'inviscid (viscMethod 0): no diffusion at all')
    # (b) モデルが構造的に適用できない
    if gas < 0 or thermal != 2 or n_species < 2:
        return out('inactive-b', 'CPG carrier' if vapor_mf > 0.0 else 'not a TP carrier (pure condensible)')
    if equil == 2:
        return out('inactive-b', 'EOS-constrained equilibrium condensation (condEquilibrium 2)')
    # (c) 実装が未対応
    if unsteady == 1 and dual == 1:
        return out('unsupported-c', 'dual-time (unsteady 1, dualTime 1)')
    if unsteady != 0 or ti != 11:
        return out('unsupported-c', 'not steady implicit pseudo-time (unsteady 0, timeIntegration 11)')
    if coupling == 2:
        return out('unsupported-c', 'speciesImplicitCoupling 2')
    if scheme != 1:
        return out('unsupported-c', 'passiveScalarScheme 0')
    if equil != 0:
        return out('unsupported-c', 'condEquilibrium 1 (relaxation form)')
    if lim_mode != 1:
        return out('unsupported-c', 'condLimiterMode 0')
    if n_cond != 1:
        return out('unsupported-c', 'nCondSpecies >= 2')
    if axisym != 0:
        return out('unsupported-c', 'axisymmetric (not verified yet)')
    if periodic:
        return out('unsupported-c', 'periodic boundary (not verified yet)')
    return out('active', 'inside the verified envelope')


def load_run(path):
    """RUN_DIR か solverConfig.yaml のパスから (y, bcond) を読む。bcondConfig.yaml が無ければ bcond は None。"""
    import yaml
    p = os.path.join(path, 'solverConfig.yaml') if os.path.isdir(path) else path
    with open(p) as f:
        y = yaml.safe_load(f) or {}
    bp = os.path.join(os.path.dirname(os.path.abspath(p)), 'bcondConfig.yaml')
    bcond = None
    if os.path.exists(bp):
        try:
            with open(bp) as f:
                bcond = yaml.safe_load(f) or {}
        except Exception:
            bcond = None
    return y, bcond


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    if not args:
        print(__doc__)
        return 2
    for a in args:
        y, b = load_run(a)
        r = classify(y, b)
        print(f'{a}: requested {r["requested"]} (default {DEFAULT}), state {r["state"]}, effective {r["effective"]} ({r["reason"]})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
